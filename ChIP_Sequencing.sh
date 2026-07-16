#!/usr/bin/env bash
# =============================================================================
# ChIP_Sequencing.sh
# -----------------------------------------------------------------------------
# Full ChIP-seq + CUT&RUN analysis pipeline for H2AK119ub1 characterisation
# in sperm and early embryos (Xenopus laevis and mouse models).
#
# Languages used:
#   Bash        — orchestration (this file)
#   Python 3    — dual-hybrid mask merging (MaskMerging.py)
#   Java        — PCR duplicate marking (Picard MarkDuplicates)
#   R           — downstream ChIP-seq statistics (ProcessingChIPdata.R)
#
# External scripts called by this pipeline:
#   MaskMerging.py          — merge PWK_PhJ x DBA_2J / C57BL_6NJ N-masked
#                             genomes into a consensus for SNPsplit
#   ProcessingChIPdata.R    — proportion tests, TF motif heatmaps,
#                             genomic characterisation
#   RECOGNICER.sh           — broad peak calling on XL9.2 (Xenopus)
#   RECOGNICER_Droso.sh     — broad peak calling on DM6 (Drosophila spike-in)
#
# Workflow:
#   1.  Directory setup + BWA hybrid genome index
#   2.  Lane merging (paired-end and single-end modes)
#   3.  Adapter/quality trimming (cutadapt + FastQC)
#   4.  BWA-MEM alignment on hybrid genome (XL9.2 + DM6 + 601DNA)
#   5.  Subgenome read extraction (XL9.2 / DM6 / 601DNA)
#   6.  QC: flagstat, fingerprint, fragment size, replicate correlation
#   7.  PCR duplicate marking + removal (Picard Java + samtools)
#   8.  Post-dedup QC
#   9.  Peak calling: MACS2 (H3K4me3) + RECOGNICER (H2Aub broad domains)
#   10. Signal tracks log2(ChIP/Input) with bamCompare
#   11. deepTools computeMatrix + plotHeatmap (Fig 2A) + plotProfile (ladder)
#   12. HOMER TF motif discovery
#   13. ChromHMM chromatin state modelling (Java)
#   14. CUT&RUN: SNPsplit genome prep + Python mask merging + bowtie2
#       + Picard dedup + SNPsplit allelic assignment
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/H2AK119ub1-sperm-embryo-exploration
# Paper   : François-Campion V. et al., Nature Communications, 2025
#           DOI: 10.1038/s41467-025-58615-7
# Date    : 2019-2022 (PhD, CR2TI UMR 1064, Nantes Universite)
# =============================================================================

#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -q max-1m.q
#$ -e ./log/
#$ -o ./log/

set -euo pipefail

# =============================================================================
# PARAMETERS — edit here
# =============================================================================

HYBRID_REF="XL9-2_dm6-40_601DNA.fa"
HYBRID_PREFIX="XL-DM-601DNA_Hybrid"
GENOME_XL="XL9_2.fa"
CHROM_SIZES="chromSizesXL9.bed"
PICARD_JAR="picard.jar"
CHROMHMM_JAR="ChromHMM.jar"

THREADS=8
DEDUP_THREADS=6
MAPQ=20

CUTADAPT_Q=20
CUTADAPT_MINLEN=20

GSIZE_XL="2.7e9"
GSIZE_DM="1.2e8"
RECOGNICER_FDR=0.001

XL_CHRS="chr1L chr1S chr2L chr2S chr3L chr3S chr4L chr4S \
          chr5L chr5S chr6L chr6S chr7L chr7S chr8L chr8S \
          chr9_10L chr9_10S"
DM_CHRS="chr_2L_fly chr_2R_fly chr_3L_fly chr_3R_fly \
          chr_4_fly chr_X_fly chr_Y_fly"


# =============================================================================
# I. DIRECTORY SETUP & BWA HYBRID GENOME INDEX
# =============================================================================

echo "[1/14] Setup and indexing..."
mkdir -p rawfastq fastq cutadapt bwa macs2 bigwig log \
         FingerPrint FragmentSize Correlation \
         HMMModel HMMClassify_All_Sample \
         bowtie SNPSplit PWK mask fasta

conda activate ChIP_SeqPipeline
bwa index -a bwtsw "${HYBRID_REF}" -p "${HYBRID_PREFIX}"


# =============================================================================
# II. LANE MERGING
# =============================================================================

echo "[2/14] Merging sequencing lanes..."

## ---- PAIRED-END ----
for sample in $(ls -1 rawfastq/*_L001_R1_001.fastq.gz \
               | cut -d "/" -f 2 | sed 's/_L001_R1_001.fastq.gz//'); do
    cat "rawfastq/${sample}_L001_R1_001.fastq.gz" \
        "rawfastq/${sample}_L002_R1_001.fastq.gz" \
        > "fastq/${sample}_1.fq.gz"
    cat "rawfastq/${sample}_L001_R2_001.fastq.gz" \
        "rawfastq/${sample}_L002_R2_001.fastq.gz" \
        > "fastq/${sample}_2.fq.gz"
done

## ---- SINGLE-END (comment out PE block above and uncomment this) ----
# for sample in $(ls -1 rawfastq/*_L001_R1_001.fastq.gz \
#               | cut -d "/" -f 2 | sed 's/_L001_R1_001.fastq.gz//'); do
#     cat "rawfastq/${sample}_L001_R1_001.fastq.gz" \
#         "rawfastq/${sample}_L002_R1_001.fastq.gz" \
#         > "fastq/${sample}_1.fq.gz"
# done


# =============================================================================
# III. ADAPTER/QUALITY TRIMMING (cutadapt + FastQC)
# =============================================================================

echo "[3/14] Trimming reads (cutadapt + FastQC)..."

## ---- PAIRED-END ----
for sample in $(ls -1 fastq/*_1.fq.gz \
               | cut -d "/" -f 2 | sed 's/_1.fq.gz//'); do
    cutadapt -q "${CUTADAPT_Q}" -m "${CUTADAPT_MINLEN}" --pair-filter=any \
        -o "cutadapt/cutadapt_${sample}_1.fastq.gz" \
        -p "cutadapt/cutadapt_${sample}_2.fastq.gz" \
        "fastq/${sample}_1.fq.gz" "fastq/${sample}_2.fq.gz"
    fastqc "cutadapt/cutadapt_${sample}_1.fastq.gz" \
           "cutadapt/cutadapt_${sample}_2.fastq.gz"
done

## ---- SINGLE-END (uncomment if needed) ----
# for sample in $(ls -1 fastq/*_1.fastq.gz \
#               | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'); do
#     cutadapt -q "${CUTADAPT_Q}" -m "${CUTADAPT_MINLEN}" \
#         -o "cutadapt/cutadapt_${sample}_1.fastq.gz" \
#         "fastq/${sample}_1.fastq.gz"
#     fastqc "cutadapt/cutadapt_${sample}_1.fastq.gz"
# done


# =============================================================================
# IV. BWA ALIGNMENT ON HYBRID GENOME (XL9.2 + DM6 + 601DNA spike-in)
# =============================================================================

echo "[4/14] BWA-MEM alignment on hybrid genome..."

## ---- PAIRED-END ----
for sample in $(ls -1 cutadapt/*_1.fastq.gz \
               | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'); do
    bwa mem -M -t "${THREADS}" "${HYBRID_PREFIX}" \
        "cutadapt/${sample}_1.fastq.gz" \
        "cutadapt/${sample}_2.fastq.gz" \
    | samtools view -hbS \
    | samtools sort \
    > "bwa/${sample}.sort.bam"
done

## ---- SINGLE-END (uncomment if needed) ----
# for sample in $(ls -1 cutadapt/*_1.fastq.gz \
#               | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'); do
#     bwa mem -M -t "${THREADS}" "${HYBRID_PREFIX}" \
#         "cutadapt/${sample}_1.fastq.gz" \
#     | samtools view -hbS | samtools sort \
#     > "bwa/${sample}.sort.bam"
# done


# =============================================================================
# V. SUBGENOME READ EXTRACTION (XL9.2 / DM6 / 601DNA)
# =============================================================================

echo "[5/14] Extracting reads per subgenome..."

for bam in bwa/*.sort.bam; do
    sample=$(basename "${bam}" .sort.bam)
    samtools index "${bam}"
    # Xenopus laevis reads
    samtools view -bh -F 0x4 "${bam}" ${XL_CHRS} \
    | samtools sort > "bwa/${sample}.XL92.sort.bam"
    # Drosophila reads (spike-in normalisation reference)
    samtools view -bh -F 0x4 "${bam}" ${DM_CHRS} \
    | samtools sort > "bwa/${sample}.DM6.sort.bam"
    # 601 DNA ladder spike-in
    samtools view -bh -F 0x4 "${bam}" 601_DNA \
    | samtools sort > "bwa/${sample}.DNA.sort.bam"
    samtools index "bwa/${sample}.XL92.sort.bam"
    samtools index "bwa/${sample}.DM6.sort.bam"
done


# =============================================================================
# VI. QC — FLAGSTAT, FINGERPRINT, FRAGMENT SIZE, CORRELATION
# =============================================================================

echo "[6/14] QC: flagstat, fingerprint, fragment size, correlation..."

# Flagstat on all BAMs
for bam in bwa/*.sort.bam; do
    sample=$(basename "${bam}" .bam)
    samtools flagstat "${bam}" > "log/${sample}.txt"
done

# Fingerprint + fragment size per sample (paired-end only)
for input_bam in bwa/*_INPUT_Input.sort.bam; do
    sample=$(basename "${input_bam}" _INPUT_Input.sort.bam)
    for suffix in ".XL92" ".DM6" ""; do
        plotFingerprint --binSize 50 -p 4 \
            -b "bwa/${sample}_H2Aub_ChIP${suffix}.sort.bam" \
               "bwa/${sample}_H3K4me3_ChIP${suffix}.sort.bam" \
               "bwa/${sample}_INPUT_Input${suffix}.sort.bam" \
            -plot "FingerPrint/${sample}_fingerprint${suffix}.svg"
        bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 \
            -b "bwa/${sample}_H2Aub_ChIP${suffix}.sort.bam" \
               "bwa/${sample}_H3K4me3_ChIP${suffix}.sort.bam" \
               "bwa/${sample}_INPUT_Input${suffix}.sort.bam" \
            -o "FragmentSize/${sample}_FragmentSize${suffix}.svg"
    done
done

# Replicate correlation (spearman)
multiBamSummary bins -p 4 --verbose \
    --bamfiles bwa/*.XL92.sort.bam \
    -o Correlation/ChipCorrelation.npz
plotCorrelation --removeOutliers \
    --corData Correlation/ChipCorrelation.npz \
    -c spearman -p heatmap --colorMap bwr --plotNumbers \
    -o Correlation/HeatmapCorChIP_spearman.svg


# =============================================================================
# VII. PCR DUPLICATE MARKING AND REMOVAL (Java — Picard MarkDuplicates)
# =============================================================================

echo "[7/14] PCR duplicate removal (Picard Java + samtools)..."

# Helper function: mark duplicates with Picard, then filter
mark_and_filter() {
    local in_bam="$1"
    local base="${in_bam%.sort.bam}"
    java -jar "${PICARD_JAR}" MarkDuplicates \
        I="${in_bam}" \
        O="${base}.duplicates.bam" \
        M="${base}.dup_metrics.txt" \
        QUIET=true
    samtools index "${base}.duplicates.bam"
    samtools view -hb -F 0x400 -q "${MAPQ}" "${base}.duplicates.bam" \
    | samtools sort -@ "${DEDUP_THREADS}" \
    > "${base}.sort.filter.bam"
    samtools index "${base}.sort.filter.bam"
}

for input_bam in bwa/*_INPUT_Input.sort.bam; do
    sample=$(basename "${input_bam}" _INPUT_Input.sort.bam)
    for genome in XL92 DM6; do
        mark_and_filter "bwa/${sample}_H3K4me3_ChIP.${genome}.sort.bam"
        mark_and_filter "bwa/${sample}_H2Aub_ChIP.${genome}.sort.bam"
        mark_and_filter "bwa/${sample}_INPUT_Input.${genome}.sort.bam"
    done
done


# =============================================================================
# VIII. POST-DEDUP QC
# =============================================================================

echo "[8/14] Post-dedup flagstat..."
for bam in bwa/*.sort.filter.bam; do
    sample=$(basename "${bam}" .bam)
    samtools flagstat "${bam}" > "log/${sample}.post_dedup.txt"
done


# =============================================================================
# IX. PEAK CALLING — MACS2 (H3K4me3 narrow) + RECOGNICER (H2Aub broad)
# NOTE: RECOGNICER requires Python 2 (separate conda env: py2)
# =============================================================================

echo "[9/14] Peak calling: MACS2 (narrow) + RECOGNICER (broad)..."

for input_bam in bwa/*_INPUT_Input.XL92.sort.filter.bam; do
    sample=$(basename "${input_bam}" _INPUT_Input.XL92.sort.filter.bam)
    mkdir -p "macs2/${sample}"

    # MACS2 — H3K4me3 narrow peaks
    for genome in XL92 DM6; do
        gsize=$([ "${genome}" = "XL92" ] && echo "${GSIZE_XL}" || echo "${GSIZE_DM}")
        macs2 callpeak \
            -t "bwa/${sample}_H3K4me3_ChIP.${genome}.sort.filter.bam" \
            -c "bwa/${sample}_INPUT_Input.${genome}.sort.filter.bam" \
            -f BAMPE --gsize="${gsize}" \
            -n "macs2/${sample}/${sample}.${genome}" \
            -B -q 0.01 --nomodel --extsize 147
    done

    # BED conversion required by RECOGNICER
    for genome in XL92 DM6; do
        bamToBed -i "bwa/${sample}_H2Aub_ChIP.${genome}.sort.filter.bam" \
            > "macs2/${sample}_H2Aub_ChIP.${genome}.sort.filter.bed"
        bamToBed -i "bwa/${sample}_INPUT_Input.${genome}.sort.filter.bam" \
            > "macs2/${sample}_INPUT_Input.${genome}.sort.filter.bed"
    done

    # RECOGNICER — H2Aub broad domains (Python 2 required)
    conda activate py2
    cd "macs2/${sample}"
    sh ../../RECOGNICER_Droso.sh \
        "../${sample}_H2Aub_ChIP.DM6.sort.filter.bed" \
        "../${sample}_INPUT_Input.DM6.sort.filter.bed" \
        "${RECOGNICER_FDR}"
    sh ../../RECOGNICER.sh \
        "../${sample}_H2Aub_ChIP.XL92.sort.filter.bed" \
        "../${sample}_INPUT_Input.XL92.sort.filter.bed" \
        "${RECOGNICER_FDR}"
    cd ../../
    conda activate ChIP_SeqPipeline

    # log2(ChIP/Input) signal tracks
    for genome in XL92 DM6; do
        suffix=$([ "${genome}" = "DM6" ] && echo "_DM6" || echo "")
        bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 \
            -b1 "bwa/${sample}_H2Aub_ChIP.${genome}.sort.filter.bam" \
            -b2 "bwa/${sample}_INPUT_Input.${genome}.sort.filter.bam" \
            -o "bigwig/${sample}_H2Aub_log2ratio${suffix}.bw"
        bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 \
            -b1 "bwa/${sample}_H3K4me3_ChIP.${genome}.sort.filter.bam" \
            -b2 "bwa/${sample}_INPUT_Input.${genome}.sort.filter.bam" \
            -o "bigwig/${sample}_H3K4me3_log2ratio${suffix}.bw"
    done
done


# =============================================================================
# X. DEEPTOOLS — computeMatrix + plotHeatmap + plotProfile
# Reproduces main figures (Fig 2A and ladder-normalised profiles)
# =============================================================================

echo "[10/14] deepTools: computeMatrix + plotHeatmap + plotProfile..."

# Fig 2A — H2Aub + H3K4me3 dynamics at USP21-sensitive gene clusters
computeMatrix reference-point -p 4 \
    --referencePoint center --verbose \
    -S  cutadapt_oo-WT_H2Aub_log2ratio.bw \
        cutadapt_oo-U21_H2Aub_log2ratio.bw \
        cutadapt_Egg-oo-WT_H2Aub_log2ratio.bw \
        cutadapt_Egg-oo-U21_H2Aub_log2ratio.bw \
        cutadapt_oo-WT_H3K4me3_log2ratio.bw \
        cutadapt_oo-U21_H3K4me3_log2ratio.bw \
        Merged_cutadapt_Egg-oo-WT_H3K4me3_log2ratio.bw \
        Merged_cutadapt_Egg-oo-U21_H3K4me3_log2ratio.bw \
    -R  Egg-WT_only.bed \
        Egg-commonU21-WT_only.bed \
        Egg-U21_only.bed \
        Egg_WT_K4Only.bed \
        Egg_commonK4_WT-U21.bed \
        Egg_U21_K4Only.bed \
    --binSize 50 --missingDataAsZero \
    --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
    --skipZeros \
    -out DeeptoolsMatrixDynamicsReplication_U21_WT.mat.gz \
    --outFileNameMatrix    DeeptoolsMatrixDynamicsReplication_U21_WT.tab \
    --outFileSortedRegions DeeptoolsMatrixDynamicsReplication_U21_WT.bed

plotHeatmap \
    --samplesLabel WTH2Aub USP21H2Aub RepWTH2Aub RepUSP21H2Aub \
                   WTK4 USP21K4 RepWTK4 RepUSP21K4 \
    --regionsLabel WT common U21 WT Common U21 \
    -m DeeptoolsMatrixDynamicsReplication_U21_WT.mat.gz \
    -out DeeptoolsMatrixDynamicsReplication_U21_WT.svg \
    --colorMap bwr \
    --outFileSortedRegions DeeptoolsMatrixDynamicsReplication_U21_WT_sorted.bed

# Ladder-normalised profiles at TSS categories
computeMatrix reference-point -p 4 \
    --referencePoint center --verbose \
    -S  cutadapt_oo-WT_LadderCov.bw \
        cutadapt_Unt_LadderCov.bw \
        cutadapt_oo-U21_LadderCov.bw \
        cutadapt_Untreated_Sp_LadderCov.bw \
        cutadapt_Egg-oo-WT_LadderCov.bw \
        cutadapt_Egg-Unt_LadderCov.bw \
        cutadapt_Egg-oo-U21_LadderCov.bw \
        cutadapt_EggExtract_Sp_LadderCov.bw \
    -R  USP21sensitivesTSS.bed MaternalOnlyTSS.bed ZygoticTSS.bed \
        MZTSS.bed OthersTSS.bed GenesTSS.bed \
    --binSize 50 --missingDataAsZero \
    --beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
    --skipZeros \
    -out Droso_USP21Sensitive_LadderH2Aub.mat.gz \
    --outFileNameMatrix    Droso_USP21Sensitive_LadderH2Aub.tab \
    --outFileSortedRegions Droso_USP21Sensitive_LadderH2Aub.region.bed

plotProfile \
    --samplesLabel \
        oo-WTH2Aub oo-UntH2Aub oo-USP21H2Aub oldo-UntH2Aub \
        egg-WTH2Aub egg-UntH2Aub egg-USP21H2Aub oldegg-UntH2Aub \
    --regionsLabel \
        USP21_DEGenes MaternalNonZygotic Zygotic MZ Others AllGenes \
    -m Droso_USP21Sensitive_LadderH2Aub.mat.gz \
    -out Droso_USP21Sensitive_LadderH2Aub.profile.svg \
    --outFileNameData Droso_USP21Sensitive_LadderH2Aub.tsv


# =============================================================================
# XI. HOMER TF MOTIF DISCOVERY
# =============================================================================

echo "[11/14] HOMER TF motif discovery..."

bedtools getfasta -fi "${GENOME_XL}" -bed Cluster3_USP21.sort.bed -tab \
    > MEME_sequences_ClusterH2Aub.bed
sed 's/chr/>chr/g' MEME_sequences_ClusterH2Aub.bed \
    | sed 's/\t/\n/g' > Sequences_H2Aub.fasta

bedtools getfasta -fi "${GENOME_XL}" -bed Background_Cluster3.bed -tab \
    > MEME_sequences_Background_ClusterH2Aub.bed
sed 's/chr/>chr/g' MEME_sequences_Background_ClusterH2Aub.bed \
    | sed 's/\t/\n/g' > MEME_sequences_Background_ClusterH2Aub.fasta

findMotifs.pl Sequences_H2Aub.fasta fasta Cluster_H2Aub \
    -fasta MEME_sequences_Background_ClusterH2Aub.fasta \
    > H2AubCluster.txt


# =============================================================================
# XII. CHROMHMM CHROMATIN STATE MODELLING (Java)
# Reference: https://github.com/jernst98/ChromHMM
# =============================================================================

echo "[12/14] ChromHMM chromatin state modelling (Java)..."
# SampleTable: BAM ChIP + Input for Spermatid, Sperm, ReplicatedSperm, Stage12

java -mx4000M -jar "${CHROMHMM_JAR}" BinarizeBam \
    -b 150 -f 2 "${CHROM_SIZES}" BAM sampleTable.tsv HMMModel

for n_states in $(seq 1 20); do
    java -mx4000M -jar "${CHROMHMM_JAR}" LearnModel \
        -p 0 -init random -b 150 -r 5000 \
        HMMModel HMMClassify_All_Sample "${n_states}" xenLae2
done


# =============================================================================
# XIII. CUT&RUN — SNPsplit genome preparation + Python mask merging
# =============================================================================

echo "[13/14] CUT&RUN: SNPsplit genome prep + Python mask merging..."
conda activate CutRun_ChenPipeline

# SNPsplit N-masked genome preparation for two strain crosses
SNPsplit_genome_preparation \
    --vcf_file SNPSplit/mgp_REL2021_snps.vcf.gz \
    --reference_genome SNPSplit \
    --dual_hybrid --strain PWK_PhJ --strain2 DBA_2J
SNPsplit_genome_preparation \
    --vcf_file SNPSplit/mgp_REL2021_snps.vcf.gz \
    --reference_genome SNPSplit \
    --dual_hybrid --strain PWK_PhJ --strain2 C57BL_6NJ

# Python 3: merge the two N-masked genomes into a single consensus FASTA
# (see MaskMerging.py for algorithm details)
# Output: mask/chrN.N-masked.fa used for bowtie2 index below
echo "  Merging N-masked genomes with MaskMerging.py (Python 3)..."
python3 MaskMerging.py \
    --dir1 PWK_PhJ_DBA_2J_dual_hybrid.based_on_GRCm39_N-masked \
    --dir2 PWK_PhJ_C57BL_6NJ_dual_hybrid.based_on_GRCm39_N-masked \
    --outdir mask

# Reformat merged FASTA (100 bp line wrapping)
for fa in mask/*.fa; do
    chrom=$(basename "${fa}" .fa)
    awk '{gsub(/.{100}/,"&\n")}1' "${fa}" > "fasta/${chrom}.fasta"
done

# Combine per-chromosome FASTA into genome file for bowtie2
cat fasta/chr*.fasta > mm39-Chen_Masked.fa

# Bowtie2 index (run once)
bowtie2-build --threads 4 -f mm39.fa mm39_gen
bowtie2-build --threads 4 -f mm39-Chen_Masked.fa mm39_genChen_masked


# =============================================================================
# XIV. CUT&RUN — Bowtie2 alignment + dedup + SNPsplit allelic assignment
# =============================================================================

echo "[14/14] CUT&RUN: Bowtie2 + Picard (Java) + SNPsplit..."

# Trimming (TrimGalore)
for sample in $(ls -1 fastq/*_1.fastq.gz \
               | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'); do
    trim_galore --fastqc --quality 20 --length 20 --paired \
        -o trimgalore \
        "fastq/${sample}_1.fastq.gz" "fastq/${sample}_2.fastq.gz"
done

# Bowtie2 alignment (CUT&RUN parameters: Chen 2021)
for sample in $(ls -1 trimgalore/*_1_val_1.fq.gz \
               | cut -d "/" -f 2 | sed 's/_1_val_1.fq.gz//'); do
    bowtie2 -p 4 --no-unal --no-mixed --no-discordant \
        -I 10 -X 700 \
        -x mm39_genChen_masked \
        -1 "trimgalore/${sample}_1_val_1.fq.gz" \
        -2 "trimgalore/${sample}_2_val_2.fq.gz" \
    | samtools view -hbS \
    | samtools sort \
    > "bowtie/${sample}.sort.bam"
    samtools index "bowtie/${sample}.sort.bam"
    samtools flagstat "bowtie/${sample}.sort.bam" > "log/${sample}.txt"
done

# Fragment size QC
for bam in bowtie/*.sort.bam; do
    sample=$(basename "${bam}" .sort.bam)
    bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 \
        -b "${bam}" -o "FragmentSize/${sample}_FragmentSize.svg"
done

# Replicate correlation (spearman)
multiBamSummary bins -bs 10000 -p 4 --verbose \
    --bamfiles \
        bowtie/H2AubEmb_Rep1.sort.bam bowtie/H2AubEmb_Rep2.sort.bam \
        bowtie/H2Aub2Cells_Rep1.sort.bam bowtie/H2Aub2Cells_Rep2.sort.bam \
        bowtie/H2AubSp_Rep1.sort.bam bowtie/H2AubSp_Rep2.sort.bam \
        bowtie/H2AubEgg_Rep1.sort.bam bowtie/H2AubEgg_Rep2.sort.bam \
    -o Correlation/CutRun_ChipCorrelation.npz
plotCorrelation --removeOutliers \
    --corData Correlation/CutRun_ChipCorrelation.npz \
    -c spearman -p heatmap --colorMap bwr --plotNumbers \
    --labels H2AubEmb1 H2AubEmb2 H2Aub2Cells1 H2Aub2Cells2 \
             H2AubSp1 H2AubSp2 H2AubOo1 H2AubOo2 \
    -o Correlation/CutRun_HeatmapCorrelation_spearman.svg

# PCR duplicate removal (Java — Picard MarkDuplicates)
for bam in bowtie/*.sort.bam; do
    sample=$(basename "${bam}" .sort.bam)
    java -jar "${PICARD_JAR}" MarkDuplicates \
        I="${bam}" \
        O="bowtie/${sample}.duplicates.bam" \
        M="bowtie/${sample}.dup_metrics.txt" \
        QUIET=true
    samtools index "bowtie/${sample}.duplicates.bam"
    samtools view -hb -F 0x400 -q "${MAPQ}" "bowtie/${sample}.duplicates.bam" \
    | samtools sort -@ "${DEDUP_THREADS}" \
    > "bowtie/${sample}.sort.filter.bam"
    samtools index "bowtie/${sample}.sort.filter.bam"
    samtools flagstat "bowtie/${sample}.sort.filter.bam" \
        > "log/${sample}.post_dedup.txt"
done

# SNPsplit allelic assignment
# Genome 1 = maternal (DBA_2J / C57BL_6NJ), Genome 2 = paternal (PWK_PhJ)
for bam in bowtie/*.sort.filter.bam; do
    sample=$(basename "${bam}" .sort.filter.bam)
    SNPsplit --paired \
        --snp_file SNPSplit/all_SNPs_PWK_PhJ_GRCm39.txt.gz \
        -o PWK "${bam}"
done

# Downstream R statistics and genomic characterisation
# (proportion tests, TF motif heatmaps — see ProcessingChIPdata.R)
# echo "  Running downstream R analysis (ProcessingChIPdata.R)..."
# Rscript ProcessingChIPdata.R

conda deactivate
echo "ChIP-seq + CUT&RUN pipeline complete."

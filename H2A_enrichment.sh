#!/usr/bin/env bash
# =============================================================================
# H2A_enrichment.sh
# -----------------------------------------------------------------------------
# H2AK119ub1 particle enrichment detection from sperm Input ChIP-seq data.
# Classifies reads by fragment size into nucleosomal / sub-nucleosomal /
# H2A-positive / H2A-negative fractions, computes per-bin coverage,
# runs proportion tests (via ProcessingChIPdata.R — standalone R script),
# discretises p-values into BED intervals, and generates spike-in normalised
# BigWig tracks and deepTools profiles.
#
# Languages used:
#   Bash    — orchestration (this file)
#   R       — proportion tests + HMM (ProcessingChIPdata.R, run separately)
#   Java    — bedGraphToBigWig is a C binary; UCSC tool
#
# Workflow:
#   1.  Fragment size classification (awk — nucleosome / H2A model)
#   2.  Per-bin coverage (bedtools intersect, 250 bp bin / 50 bp slide)
#   3.  NOTE: proportion tests run via: Rscript ProcessingChIPdata.R
#   4.  P-value discretisation + BED merging at multiple FDR thresholds
#   5.  Per-chromosome coverage of enriched/depleted regions
#   6.  Intersection with H2AK119ub1 peaks
#   7.  Fragment subset extraction + BAM generation per fraction
#   8.  Fragment size QC per fraction (bamPEFragmentSize)
#   9.  log2 ratio BigWig tracks (bamCompare)
#   10. 601 DNA ladder spike-in normalisation + BigWig
#   11. deepTools computeMatrix + plotProfile (ladder-normalised profiles)
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/H2AK119ub1-sperm-embryo-exploration
# Paper   : François-Campion V. et al., Nature Communications, 2025
#           DOI: 10.1038/s41467-025-58615-7
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

SAMPLE="Input_Sp_WithoutDupl"
GENOME_BIN="Genome/XLaevis_250bpbin_slid50.bed"
CHROM_BED="GFF/XL9_Chromosomes.bed"
CHROM_SIZES="chromSizesXL9.bed"
PEAKS="Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed"

# 601 DNA spike-in counts — update from idxstats output
CHIP_601_READS=10807
INPUT_601_READS=484
LADDER_RATIO=$(echo "scale=6; ${INPUT_601_READS}/${CHIP_601_READS}" | bc)

THREADS=4
MAPQ=20
MIN_FRAG_READS=5     # min total reads per bin for analysis

# FDR thresholds for BED discretisation
PVAL_LABELS=(0.1 0.05 0.001)

# =============================================================================
# I. DIRECTORY SETUP
# =============================================================================

mkdir -p Merged Fragment IGV Coverage DNALadder FragmentSize

# =============================================================================
# II. FRAGMENT SIZE CLASSIFICATION (awk)
# =============================================================================

echo "[1/11] Classifying reads by fragment size (TLEN)..."

# Convert sorted BAM to BED with absolute fragment length in col 5
samtools view "Merged/${SAMPLE}.sort.bam" \
    | awk '{print $3"\t"$4"\t"$4+50"\t"$1"\t"sqrt(($9)^2)}' \
    > "Merged/${SAMPLE}.sort.bed"

# Nucleosome model: >135 bp = nucleosome(1), 31-135 = sub-nucleosomal(2)
awk '{
    if      ($5 > 135)                   print $0"\t1"
    else if ($5 >  30 && $5 <= 135)      print $0"\t2"
    else                                 print $0"\t0"
}' "Merged/${SAMPLE}.sort.bed" > "Merged/${SAMPLE}_flagNucl.bed"

# H2A model: >95 bp = H2A-positive(1), 31-95 = H2A-negative(2)
awk '{
    if      ($5 > 95)                    print $0"\t1"
    else if ($5 >  30 && $5 <= 95)       print $0"\t2"
    else                                 print $0"\t0"
}' "Merged/${SAMPLE}.sort.bed" > "Merged/${SAMPLE}_H2A.bed"

# Split nucleosome model into two BED files
for flag_val in 1 2; do
    label=$([ "${flag_val}" -eq 1 ] && echo "Nucl" || echo "SubNucl")
    grep $'\t'"${flag_val}"'$' "Merged/${SAMPLE}_flagNucl.bed" \
        | sort -k1,1 -k2,2n \
        > "Merged/${SAMPLE}_${label}.bed"
done

# Split H2A model into two BED files
for flag_val in 1 2; do
    label=$([ "${flag_val}" -eq 1 ] && echo "H2APos" || echo "H2ANeg")
    grep $'\t'"${flag_val}"'$' "Merged/${SAMPLE}_H2A.bed" \
        | sort -k1,1 -k2,2n \
        > "Merged/${SAMPLE}_${label}.bed"
done


# =============================================================================
# III. PER-BIN COVERAGE (bedtools intersect — 250 bp / slide 50 bp)
# =============================================================================

echo "[2/11] Computing read coverage per genomic bin..."

for fraction in Nucl SubNucl H2APos H2ANeg; do
    bedtools intersect \
        -sorted \
        -a "${GENOME_BIN}" \
        -b "Merged/${SAMPLE}_${fraction}.bed" \
        -c \
        > "Fragment/${SAMPLE}_${fraction}_slid50.bed"
done

# Combine Nucl + SubNucl into 6-column bedgraph (col4=Nucl, col5=SubNucl, col6=Total)
paste "Fragment/${SAMPLE}_Nucl_slid50.bed" \
      "Fragment/${SAMPLE}_SubNucl_slid50.bed" \
    | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"($4+$8)}' \
    > "Fragment/${SAMPLE}.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_Nucl_vs_SemiNucl.full.bedgraph"

# Combine H2APos + H2ANeg
paste "Fragment/${SAMPLE}_H2APos_slid50.bed" \
      "Fragment/${SAMPLE}_H2ANeg_slid50.bed" \
    | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"($4+$8)}' \
    > "Fragment/${SAMPLE}.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_H2A_vs_NonH2A.full.bedgraph"

# Filter bins with >= MIN_FRAG_READS total reads
for bg in \
    "Fragment/${SAMPLE}.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_Nucl_vs_SemiNucl.full.bedgraph" \
    "Fragment/${SAMPLE}.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_H2A_vs_NonH2A.full.bedgraph"; do
    awk -v min="${MIN_FRAG_READS}" '$6 >= min' "${bg}" \
        > "${bg%.full.bedgraph}_above5.full.bedgraph"
done


# =============================================================================
# NOTE: PROPORTION TESTS — run separately via R
# =============================================================================
# The .probabilities output files required by Section IV below are generated
# by running:
#   Rscript ProcessingChIPdata.R
# (Section VI of that script — standalone, not called from here)


# =============================================================================
# IV. P-VALUE DISCRETISATION AND BED MERGING
# =============================================================================

echo "[3/11] Discretising p-values and merging BED intervals..."

H2A_PROB="Fragment/${SAMPLE}.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities"
NUCL_PROB="Fragment/${SAMPLE}.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities"

# Discretise H2A model: 3 levels (p<=0.1=1, p<=0.05=2, p<=0.001=3)
awk '{
    p1=0; p2=0
    if ($6 > 0) {
        if ($9  <= 0.1)   p1=1
        if ($9  <= 0.05)  p1=2
        if ($9  <= 0.001) p1=3
        if ($10 <= 0.1)   p2=1
        if ($10 <= 0.05)  p2=2
        if ($10 <= 0.001) p2=3
    }
    print $0"\t"p1"\t"p2
}' "${H2A_PROB}" > "${H2A_PROB}.discrete"

# Discretise Nucl model: 5 levels (adds p<=0.0001 and p<=0.00001)
awk '{
    p1=0; p2=0
    if ($6 > 0) {
        if ($9  <= 0.1)      p1=1
        if ($9  <= 0.05)     p1=2
        if ($9  <= 0.001)    p1=3
        if ($9  <= 0.0001)   p1=4
        if ($9  <= 0.00001)  p1=5
        if ($10 <= 0.1)      p2=1
        if ($10 <= 0.05)     p2=2
        if ($10 <= 0.001)    p2=3
        if ($10 <= 0.0001)   p2=4
        if ($10 <= 0.00001)  p2=5
    }
    print $0"\t"p1"\t"p2
}' "${NUCL_PROB}" > "${NUCL_PROB}.discrete"

# =============================================================================
# NOTE: HIDDEN MARKOV MODEL — run separately via R
# =============================================================================
# The _hmm_state_50bpw.full output files not required by Section IV below are generated
# by running:
#   Rscript ProcessingChIPdata.

# =============================================================================
# IV. P-VALUE DISCRETISATION AND BED MERGING
# =============================================================================

# Merge BED at each FDR threshold
for label in "${PVAL_LABELS[@]}"; do
    lvl=$(echo "${label}" | awk '{if($1==0.1) print 1; else if($1==0.05) print 2; else print 3}')
    awk -v t="${lvl}" '$11 >= t' "${H2A_PROB}.discrete"  | mergeBed > "IGV/H2Aenriched_${label}.bed"
    awk -v t="${lvl}" '$12 >= t' "${H2A_PROB}.discrete"  | mergeBed > "IGV/H2Adepleted_${label}.bed"
    awk -v t="${lvl}" '$11 >= t' "${NUCL_PROB}.discrete" | mergeBed > "IGV/NucleosomeEnriched_${label}.bed"
    awk -v t="${lvl}" '$12 >= t' "${NUCL_PROB}.discrete" | mergeBed > "IGV/NucleosomeDepleted_${label}.bed"
done


# =============================================================================
# V. PER-CHROMOSOME COVERAGE OF ENRICHED/DEPLETED REGIONS
# =============================================================================

echo "[4/11] Computing per-chromosome coverage..."

for bed in IGV/*.bed; do
    name=$(basename "${bed}" .bed)
    bedtools coverage \
        -a "${CHROM_BED}" \
        -b "${bed}" \
        > "Coverage/Coverage_${name}.bed"
done


# =============================================================================
# VI. INTERSECTION WITH H2AK119ub1 PEAKS
# =============================================================================

echo "[5/11] Intersecting with H2AK119ub1 peaks..."

for label in "${PVAL_LABELS[@]}"; do
    for type in NucleosomeEnriched NucleosomeDepleted H2Aenriched H2Adepleted; do
        bedtools coverage \
            -a "${PEAKS}" \
            -b "IGV/${type}_${label}.bed" \
            > "Coverage/${type}_${label}.coverage"
    done
done


# =============================================================================
# VII. FRAGMENT SUBSET EXTRACTION — BAM per fraction
# =============================================================================

echo "[6/11] Extracting BAM per fragment fraction..."

for type in Nucl SubNucl; do
    sort -k1,1 -k2,2n "Merged/${SAMPLE}_${type}.bed" \
        | mergeBed \
        > "Merged/${SAMPLE}_${type}.sorted.bed"
done
for type in H2APos H2ANeg; do
    sort -k1,1 -k2,2n "Merged/${SAMPLE}_${type}.bed" \
        | mergeBed \
        > "Merged/${SAMPLE}_${type}.sorted.bed"
done

declare -A FRACTION_BED=(
    [Nucleosome]="Merged/${SAMPLE}_Nucl.bed"
    [SubNucleosome]="Merged/${SAMPLE}_SubNucl.sorted.bed"
    [WithH2A]="Merged/${SAMPLE}_H2APos.sorted.bed"
    [WithoutH2A]="Merged/${SAMPLE}_H2ANeg.sorted.bed"
)

for label in "${!FRACTION_BED[@]}"; do
    bedtools intersect \
        -a "Merged/${SAMPLE}.sort.bam" \
        -b "${FRACTION_BED[$label]}" \
        > "IGV/Input_${label}.bam"
    samtools sort -@ 6 "IGV/Input_${label}.bam" \
        -o "IGV/Input_${label}.sort.bam"
    samtools index "IGV/Input_${label}.sort.bam"
done


# =============================================================================
# VIII. FRAGMENT SIZE QC PER FRACTION
# =============================================================================

echo "[7/11] Fragment size QC per fraction..."

for bam in IGV/Input_*.sort.bam; do
    sample=$(basename "${bam}" .sort.bam)
    bamPEFragmentSize \
        --maxFragmentLength 300 --binSize 1000 \
        -p "${THREADS}" \
        -b "${bam}" \
        -T "Fragment size — ${sample}" \
        --plotFileFormat pdf \
        -o "FragmentSize/${sample}_FragmentSize.pdf" \
        --table "FragmentSize/${sample}.tsv" \
        --samplesLabel "${sample}"
done

# MidOverlap fraction (if pre-generated)
for bam in FragmentSize/*_MidOverlap.bam; do
    [ -f "${bam}" ] || continue
    sample=$(basename "${bam}" .bam)
    samtools sort  -@ 6 "${bam}" -o "FragmentSize/${sample}.sort.bam"
    samtools index "FragmentSize/${sample}.sort.bam"
    bamPEFragmentSize \
        -p "${THREADS}" --plotFileFormat pdf \
        --maxFragmentLength 300 \
        -b "FragmentSize/${sample}.sort.bam" \
        -T "Fragment size" \
        -o "FragmentSize/${sample}_Hist.pdf" \
        --table "FragmentSize/${sample}.tsv" \
        --samplesLabel "${sample}"
done


# =============================================================================
# IX. LOG2 RATIO BIGWIG TRACKS (bamCompare)
# =============================================================================

echo "[8/11] log2(Nucleosome/SubNucleosome) and log2(H2APos/H2ANeg) BigWig..."

bamCompare --binSize 50 -p "${THREADS}" -v -of bigwig --operation log2 \
    -b1 IGV/Input_Nucleosome.sort.bam \
    -b2 IGV/Input_SubNucleosome.sort.bam \
    -o IGV/NuclPos_vs_SubNuclNeg_log2ratio.bw

bamCompare --binSize 50 -p "${THREADS}" -v -of bigwig --operation log2 \
    -b1 IGV/Input_WithH2A.sort.bam \
    -b2 IGV/Input_WithoutH2A.sort.bam \
    -o IGV/H2APos_vs_H2ANeg_log2ratio.bw


# =============================================================================
# X. 601 DNA LADDER SPIKE-IN NORMALISATION → BigWig
# Normalisation ratio = INPUT_601_reads / CHIP_601_reads
# =============================================================================

echo "[9/11] 601 DNA ladder spike-in normalisation (ratio = ${LADDER_RATIO})..."

# 50 bp bin version (trim 100 bp from each end of 250 bp bins)
awk '{print $1"\t"($2+100)"\t"($3-100)}' \
    "DNALadder/XLaevis_250bpbin_slid50.bed" \
    > "DNALadder/XLaevis_50bpbin.bed"

for input_bam in DNALadder/*_INPUT_Input.XL92.sort.filter.bam; do
    sample=$(basename "${input_bam}" _INPUT_Input.XL92.sort.filter.bam)

    bedtools intersect \
        -a "DNALadder/XLaevis_250bpbin.bed" \
        -b "DNALadder/${sample}_INPUT_Input.XL92.sort.filter.bam" -c \
        > "DNALadder/${sample}_INPUT_Input.GenomeCov.bed"

    bedtools intersect \
        -a "DNALadder/XLaevis_250bpbin.bed" \
        -b "DNALadder/${sample}_H2Aub_ChIP.XL92.sort.filter.bam" -c \
        > "DNALadder/${sample}_H2Aub_ChIP.GenomeCov.bed"

    # Ladder normalisation: (ChIP_reads / (Input_reads * 2)) / ratio
    paste "DNALadder/${sample}_H2Aub_ChIP.GenomeCov.bed" \
          "DNALadder/${sample}_INPUT_Input.GenomeCov.bed" \
        | awk -v ratio="${LADDER_RATIO}" '{
            if ($8 > 0)
                print $1"\t"$2"\t"$3"\t"$4"\t"($8*2)"\t"($4/(($8)*2))"\t"($4/(($8)*2))/ratio
            else
                print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t0\t0"
        }' > "DNALadder/${sample}_LadderCov.bed"

    awk '{print $1"\t"$2"\t"$3"\t"$7}' "DNALadder/${sample}_LadderCov.bed" \
        > "DNALadder/${sample}_LadderCov.bg"

    ./bedGraphtoBigWig \
        "DNALadder/${sample}_LadderCov.bg" \
        "${CHROM_SIZES}" \
        "DNALadder/${sample}_LadderCov.bw"
    echo "  BigWig saved: ${sample}_LadderCov.bw"
done


# =============================================================================
# XI. DEEPTOOLS — computeMatrix + plotProfile (ladder-normalised profiles)
# =============================================================================

echo "[10/11] deepTools computeMatrix + plotProfile (ladder profiles)..."

computeMatrix reference-point \
    -p "${THREADS}" \
    --referencePoint center --verbose \
    -S  cutadapt_oo-WT_LadderCov.bw \
        cutadapt_Unt_LadderCov.bw \
        cutadapt_oo-U21_LadderCov.bw \
        cutadapt_Untreated_Sp_LadderCov.bw \
        cutadapt_Egg-oo-WT_LadderCov.bw \
        cutadapt_Egg-Unt_LadderCov.bw \
        cutadapt_Egg-oo-U21_LadderCov.bw \
        cutadapt_EggExtract_Sp_LadderCov.bw \
    -R  USP21sensitivesTSS.bed \
        MaternalOnlyTSS.bed \
        ZygoticTSS.bed \
        MZTSS.bed \
        OthersTSS.bed \
        GenesTSS.bed \
    --binSize 50 --missingDataAsZero \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength  5000 \
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

echo "[11/11] H2A enrichment pipeline complete."

#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -q max-1m.q
#$ -e  ./log/
#$ -o ./log/

bwa index -a bwtsw XL9-2_dm6-40_601DNA.fa -p XL-DM-601DNA_Hybrid
conda activate ChIP_SeqPipeline
mkdir rawfastq
mkdir fastq
mkdir FingerPrint
mkdir FragmentSize
mkdir bwa
mkdir cutadapt
mkdir macs2
mkdir log


## Put your FASTQ in rawfastq folder
## RENAME your sample like : <SAMPLE_NAME><HISTONEMARK_ChIP or INPUT_Input><SEQUENCINGLANEparameter><fastq.gz>
# Example : Sperm1_H2AK119ub1_ChIP_L001_R1_001.fastq.gz and Sperm1_INPUT_Input_L001_R1_001.fastq.gz

############################
### FOR PAIRED-END DATA ####
############################

# Merge Lanes
for f in `ls -1 rawfastq/*_L001_R1_001.fastq.gz | cut -d "/" -f 2 | sed 's/_L001_R1_001.fastq.gz//'`
	do 
	cat rawfastq/${f}_L001_R1_001.fastq.gz rawfastq/${f}_L002_R1_001.fastq.gz > fastq/${f}_1.fq.gz
	cat rawfastq/${f}_L001_R2_001.fastq.gz rawfastq/${f}_L002_R2_001.fastq.gz > fastq/${f}_2.fq.gz
done

# low quality reads trimming
for f in `ls -1 fastq/*_1.fq.gz | cut -d "/" -f 2 | sed 's/_1.fq.gz//'`
	do cutadapt -q 20 -m 20 --pair-filter=any -o cutadapt/cutadapt_${f}_1.fastq.gz -p cutadapt/cutadapt_${f}_2.fastq.gz fastq/${f}_1.fq.gz fastq/${f}_2.fq.gz
	fastqc cutadapt/cutadapt_${f}_1.fastq.gz cutadapt/cutadapt_${f}_2.fastq.gz
done

# Alignment PE
for f in `ls -1 cutadapt/*_1.fastq.gz  | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do bwa mem -M -t 8 ~/../../../mnt/c/Documents\ and\ Settings/ValentinFC/Documents/These/Analyses/GFF/XL-DM-601DNA_Hybrid cutadapt/${f}_1.fastq.gz cutadapt/${f}_2.fastq.gz | samtools view -hbS | samtools sort > bwa/${f}.sort.bam
done

############################
### FOR SINGLE-END DATA ####
############################

# Merge Lanes
for f in `ls -1 rawfastq/*_L001_R1_001.fastq.gz | cut -d "/" -f 2 | sed 's/_L001_R1_001.fastq.gz//'`
	do 
	cat rawfastq/${f}_L001_R1_001.fastq.gz rawfastq/${f}_L002_R1_001.fastq.gz > fastq/${f}_1.fq.gz
done

# Low quality reads trimming
for f in `ls -1 fastq/*_1.fastq.gz  | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do cutadapt -q 20 -m 20 -o cutadapt/cutadapt_${f}_1.fastq.gz fastq/${f}_1.fastq.gz 
	fastqc cutadapt/cutadapt_${f}_1.fastq.gz
done

##########################################
### Reads alignment into hybrid genome ###
##########################################

for f in `ls -1 cutadapt/*_1.fastq.gz  | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do bwa mem -M -t 8 ~/../../../mnt/c/Documents\ and\ Settings/ValentinFC/Documents/These/Analyses/GFF/XL-DM-601DNA_Hybrid cutadapt/${f}_1.fastq.gz | samtools view -hbS | samtools sort > bwa/${f}.sort.bam
done

###################################
### READ EXTRACTION PER GENOME ####
###################################

# Index BAM file, and extract subgenomes chromosomes (601 DNA ladder, Drosophila melanogaster, Xenopus laevis)
for f in `ls -1 bwa/*.bam  | cut -d "/" -f 2 | sed 's/.bam//'`
	do samtools index bwa/${f}.sort.bam 
	   samtools view -bh -F 0x4 bwa/${f}.sort.bam chr1L chr1S chr2L chr2S chr3L chr3S chr4L chr4S chr5L chr5S chr6L chr6S chr7L chr7S chr8L chr8S chr9_10L chr9_10S | samtools view -bS | samtools sort > bwa/${f}.XL92.sort.bam 
	   samtools view -bh -F 0x4 bwa/${f}.sort.bam chr_2L_fly chr_2R_fly chr_3L_fly chr_3R_fly chr_4_fly chr_X_fly chr_Y_fly | samtools view -bS | samtools sort > bwa/${f}.DM6.sort.bam
	   samtools view -bh -F 0x4 bwa/${f}.sort.bam 601_DNA | samtools view -bS | samtools sort > bwa/${f}.DNA.sort.bam
	   samtools index bwa/${f}.XL92.sort.bam 
	   samtools index bwa/${f}.DM6.sort.bam
done

################
### DATA QC ####
################

#Sample quality
for f in `ls -1 bwa/*sort.bam | cut -d "/" -f 2 | sed 's/.bam//'`
	do samtools flagstat bwa/${f}.bam > log/${f}.txt
done

#Replicates quality
multiBamSummary bins -p 4 --verbose --bamfiles <ALL_BAMFILES> -o Correlation/ChipCorrelationChIP_Input.npz
plotCorrelation --removeOutliers --corData Correlation/ChipCorrelationChIP_Input.npz -c spearman -p heatmap --colorMap bwr --plotNumbers --labels <ALL_SamplesLABEL> -o Correlation/HeatmapCorChIP_Inputspearman.svg

########################
### Histone mark QC ####
########################

# Fragmentsize only on PAIRED-END DATA
for f in `ls -1 bwa/*_INPUT_Input.sort.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.sort.bam//'`
	do 
	   plotFingerprint --binSize 50 -p 4 -b bwa/${f}_H2Aub_ChIP.XL92.sort.bam bwa/${f}_H3K4me3_ChIP.XL92.sort.bam bwa/${f}_INPUT_Input.XL92.sort.bam -plot FingerPrint/${f}_fingerprint.XL92.svg 
	   plotFingerprint --binSize 50 -p 4 -b bwa/${f}_H2Aub_ChIP.DM6.sort.bam bwa/${f}_H3K4me3_ChIP.DM6.sort.bam bwa/${f}_INPUT_Input.DM6.sort.bam -plot FingerPrint/${f}_fingerprint.DM6.svg
	   plotFingerprint --binSize 50 -p 4 -b bwa/${f}_H2Aub_ChIP.sort.bam bwa/${f}_H3K4me3_ChIP.sort.bam bwa/${f}_INPUT_Input.sort.bam -plot FingerPrint/${f}_fingerprint.svg
	   
	   ### ONLY ON PAIRED-END DATA
	   bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bwa/${f}_H2Aub_ChIP.XL92.sort.bam bwa/${f}_H3K4me3_ChIP.XL92.sort.bam bwa/${f}_INPUT_Input.XL92.sort.bam -o FragmentSize/${f}_FragmentSize.XL92.svg
	   bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bwa/${f}_H2Aub_ChIP.DM6.sort.bam bwa/${f}_H3K4me3_ChIP.DM6.sort.bam bwa/${f}_INPUT_Input.DM6.sort.bam -o FragmentSize/${f}_FragmentSize.DM6.svg
	   bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bwa/${f}_H2Aub_ChIP.sort.bam bwa/${f}_H3K4me3_ChIP.sort.bam bwa/${f}_INPUT_Input.sort.bam -o FragmentSize/${f}_FragmentSize.svg
done

###############################
### PCR duplicates removal ####
###############################

# mark duplicates with Picard && index result bam for DM6 and XL92 genome, on each <HISTONEMARK> and <INPUT>
for f in `ls -1 bwa/*_INPUT_Input.sort.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.sort.bam//'`
	do 
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H3K4me3_ChIP.XL92.sort.bam O=bwa/${f}_H3K4me3_ChIP.XL92.duplicates.bam M=bwa/${f}_H3K4me3_ChIP.XL92.dup_metrics.txt
	   samtools index bwa/${f}_H3K4me3_chIP.XL92.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H2Aub_ChIP.XL92.sort.bam O=bwa/${f}_H2Aub_ChIP.XL92.duplicates.bam M=bwa/${f}_H2Aub_ChIP.XL92.dup_metrics.txt
	   samtools index bwa/${f}_H2Aub_chIP.XL92.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_INPUT_Input.XL92.sort.bam O=bwa/${f}_INPUT_Input.XL92.duplicates.bam M=bwa/${f}_INPUT_Input.XL92.dup_metrics.txt
	   samtools index bwa/${f}_INPUT_Input.XL92.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H3K4me3_ChIP.DM6.sort.bam O=bwa/${f}_H3K4me3_ChIP.DM6.duplicates.bam M=bwa/${f}_H3K4me3_ChIP.DM6.dup_metrics.txt
	   samtools index bwa/${f}_H3K4me3_chIP.DM6.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_H2Aub_ChIP.DM6.sort.bam O=bwa/${f}_H2Aub_ChIP.DM6.duplicates.bam M=bwa/${f}_H2Aub_ChIP.DM6.dup_metrics.txt
	   samtools index bwa/${f}_H2Aub_chIP.DM6.duplicates.bam
	   java -jar picard.jar MarkDuplicates I=bwa/${f}_INPUT_Input.DM6.sort.bam O=bwa/${f}_INPUT_Input.DM6.duplicates.bam M=bwa/${f}_INPUT_Input.DM6.dup_metrics.txt
	   samtools index bwa/${f}_INPUT_Input.DM6.duplicates.bam
done	   

# filter out duplicates
for f in `ls -1 bwa/*.duplicates.bam | cut -d "/" -f 2 | sed 's/.duplicates.bam//'`
	do samtools view -hb -F 0x400 -q 20 bwa/${f}.duplicates.bam | samtools sort -@ 6 > bwa/${f}.sort.filter.bam
	   samtools index bwa/${f}.sort.filter.bam 
done

#############################################
### DATA QC after PCR duplicates removal ####
#############################################

for f in `ls -1 bwa/*.sort.filter.bam   | cut -d "/" -f 2 | sed 's/.bam//'`
	do samtools flagstat bwa/${f}.bam > log/${f}.txt
done

#############################################
### PEAKS calling and DATA visualisation ####
#############################################

for f in `ls -1 bwa/*_INPUT_Input.XL92.sort.filter.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.XL92.sort.filter.bam//'`
	do 
	mkdir macs2/${f}
 	# peak calling for H3K4me3
	macs2 callpeak -t bwa/${f}_H3K4me3_ChIP.DM6.sort.filter.bam -c bwa/${f}_INPUT_Input.DM6.sort.filter.bam -f BAMPE --gsize=1.2e8 -n macs2/${f}/${f}.DM6 -B -q 0.01 --nomodel --extsize 147
	macs2 callpeak -t bwa/${f}_H3K4me3_ChIP.XL92.sort.filter.bam -c bwa/${f}_INPUT_Input.XL92.sort.filter.bam -f BAMPE --gsize=2.6e9 -n macs2/${f}/${f}.XL92 -B -q 0.01 --nomodel --extsize 147

 	# setup for RECOGNICER peak calling (bedfiles)
	bamToBed -i bwa/${f}_H2Aub_ChIP.DM6.sort.filter.bam > macs2/${f}_H2Aub_chIP.DM6.sort.filter.bed
	bamToBed -i bwa/${f}_INPUT_Input.DM6.sort.filter.bam > macs2/${f}_INPUT_Input.DM6.sort.filter.bed
	bamToBed -i bwa/${f}_H2Aub_ChIP.XL92.sort.filter.bam > macs2/${f}_H2Aub_chIP.XL92.sort.filter.bed
	bamToBed -i bwa/${f}_INPUT_Input.XL92.sort.filter.bam > macs2/${f}_INPUT_Input.XL92.sort.filter.bed

 	# Signal track log(ChIP/Input) of H2AK119ub1 and H3K4me3
	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H2Aub_chIP.XL92.sort.filter.bam -b2 bwa/${f}_INPUT_Input.XL92.sort.filter.bam -o bigwig/${f}_H2Aub_log2ratio.bw
	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H3K4me3_ChIP.XL92.sort.filter.bam -b2 bwa/${f}_INPUT_Input.XL92.sort.filter.bam -o bigwig/${f}_H3K4me3_log2ratio.bw
	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H2Aub_chIP.DM6.sort.filter.bam -b2 bwa/${f}_INPUT_Input.DM6.sort.filter.bam -o bigwig/${f}_H2Aub_log2ratio_DM6.bw
	bamCompare --binSize 50 -p 4 -v -of bigwig --operation log2 -b1 bwa/${f}_H3K4me3_ChIP.DM6.sort.filter.bam -b2 bwa/${f}_INPUT_Input.DM6.sort.filter.bam -o bigwig/${f}_H3K4me3_log2ratio_DM6.bw

 	# peak calling with RECOGNICER for H2AK119ub1 peaks
	# activate py2 env for RECOGNICER
	conda activate py2
	cd macs2/${f}
	sh ../../RECOGNICER_Droso.sh ../${f}_H2Aub_chIP.DM6.sort.filter.bed ../${f}_INPUT_Input.DM6.sort.filter.bed 0.001
	sh ../../RECOGNICER.sh ../${f}_H2Aub_chIP.XL92.sort.filter.bed ../${f}_INPUT_Input.XL92.sort.filter.bed 0.001
	cd ../../
	conda activate ChIP_SeqPipeline
done

#Fig 2A example
computeMatrix reference-point -p 4 --referencePoint center --verbose -S cutadapt_oo-WT_H2Aub_log2ratio.bw cutadapt_oo-U21_H2Aub_log2ratio.bw cutadapt_Egg-oo-WT_H2Aub_log2ratio.bw cutadapt_Egg-oo-U21_H2Aub_log2ratio.bw cutadapt_oo-WT_H3K4me3_log2ratio.bw cutadapt_oo-U21_H3K4me3_log2ratio.bw Merged_cutadapt_Egg-oo-WT_H3K4me3_log2ratio.bw Merged_cutadapt_Egg-oo-U21_H3K4me3_log2ratio.bw -R Egg-WT_only.bed Egg-commonU21-WT_only.bed Egg-U21_only.bed Egg_WT_K4Only.bed Egg_commonK4_WT-U21.bed Egg_U21_K4Only.bed --binSize 50 --missingDataAsZero --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --skipZeros -out DeeptoolsMatrixDynamicsReplication_U21_WT.mat.gz --outFileNameMatrix DeeptoolsMatrixDynamicsReplication_U21_WT.tab --outFileSortedRegions DeeptoolsMatrixDynamicsReplication_U21_WT.bed
plotHeatmap --samplesLabel WTH2Aub USP21H2Aub RepWTH2Aub RepUSP21H2Aub WTK4 USP21K4 RepWTK4 RepUSP21K4 --regionsLabel WT common U21 WT Common U21 -m DeeptoolsMatrixDynamicsReplication_U21_WT.mat.gz -out DeeptoolsMatrixDynamicsReplication_U21_WT.svg --colorMap bwr --outFileSortedRegions DeeptoolsMatrixDynamicsReplication_U21_WT.bed

#######################################
### HOMER MOTIF DISCOVERY ON PEAKS ####
#######################################

# extract fasta sequences from a bedfile coordinates of interest
bedtools getfasta -fi XL9_2.fa -bed Cluster3_USP21.sort.bed -tab > MEME_sequences_ClusterH2Aub.bed
sed 's/chr/>chr/g' MEME_sequences_ClusterH2Aub.bed | sed 's/\t/\n/g' > Sequences_H2Aub.fasta

# extract fasta sequences from a background bedfile coordinates
bedtools getfasta -fi XL9_2.fa -bed Background_Cluster3.bed -tab > MEME_sequences_Background_ClusterH2Aub.bed
sed 's/chr/>chr/g' MEME_sequences_Background_ClusterH2Aub.bed | sed 's/\t/\n/g' > MEME_sequences_Background_ClusterH2Aub.fasta

# homer motif research command line
findMotifs.pl Sequences_H2Aub.fasta fasta Cluster_H2Aub -fasta MEME_sequences_Background_ClusterH2Aub.fasta > H2AubCluster.txt

#######################
### CHROMHMM MODEL ####
#######################

https://github.com/jernst98/ChromHMM/blob/master/ChromHMM_manual.pdf
mkdir HMMModel
mkdir HMMClassify

#SampleTable contains all BAM ChIP and Input samples shown in the ChromHMM model : Spermatid, Sperm, ReplicatedSperm and Stage12 H2Aub
java -mx4000M -jar ChromHMM.jar BinarizeBam -b 150 -f 2 ~/../../../mnt/c/Documents\ and\ Settings/ValentinFC/Documents/These/Analyses/GFF/chromSizesXL9.bed BAM sampleTable.tsv HMMModel
#Model generated for 0 at 20 clusters
for f in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
	do java -mx4000M -jar ChromHMM.jar LearnModel -p 0 -init random -b 150 -r 5000 HMMModel HMMClassify_All_Sample ${f} xenLae2
done

################################
#### CUT&RUN CHEN pipeline #####
################################

conda activate CutRun_ChenPipeline
#SNPSplit genome extraction : Paternal = PWK_PhJ, Mat = DBA_2JxC57BL_6NJ https://felixkrueger.github.io/SNPsplit/genome_prep/genome_preparation/
SNPsplit_genome_preparation --vcf_file SNPSplit/mgp_REL2021_snps.vcf.gz --reference_genome SNPSplit --dual_hybrid --strain PWK_PhJ --strain2 DBA_2J
SNPsplit_genome_preparation --vcf_file SNPSplit/mgp_REL2021_snps.vcf.gz --reference_genome SNPSplit --dual_hybrid --strain PWK_PhJ --strain2 C57BL_6NJ

MaskMerging.py PWK_PhJ_DBA_2J_N-masked.fa PWK_PhJ_C57BL_6NJ_N-masked.fa

for f in `ls -1 mask/*.fa | cut -d "/" -f 2 | sed 's/.fa//'`
	do awk '{gsub(/.{100}/,"&\n")}1' mask/${f}.fa > fasta/${f}.fasta 
done

cat fasta/*.fasta > mm39-Chen_Masked.fa
bowtie2-build --threads 4 -f mm39.fa mm39_gen
bowtie2-build --threads 4 -f mm39-Chen_Masked.fa mm39_genChen_masked

#Remove adapter/low quality reads and export fastqc report
for f in `ls -1 fastq/*_1.fastq.gz | cut -d "/" -f 2 | sed 's/_1.fastq.gz//'`
	do trim_galore --fastqc --quality 20 --length 20 --paired -o trimgalore fastq/${f}_1.fastq.gz fastq/${f}_2.fastq.gz
done

# bowtie 2 alignment for CUT&RUN data parameters on bibliography
for f in `ls -1 trimgalore/*_1.fq.gz  | cut -d "/" -f 2 | sed 's/_1_val_1.fq.gz//'`
	do bowtie2 -p 4 --no-unal --no-mixed --no-discordant -I 10 -X 700 -x ~/../../../mnt/c/Documents\ and\ Settings/ValentinFC/Documents/These/Analyses/GFF/mm39_genChen_masked -1 trimgalore/${f}_1_val_1.fq.gz -2 trimgalore/${f}_2_val_2.fq.gz | samtools view -hbS | samtools sort > bowtie/${f}.sort.bam
done

# Check fragmentsize, and compute alignement data
for f in `ls -1 bowtie/*.sort.bam  | cut -d "/" -f 2 | sed 's/.sort.bam//'`
	do 
	samtools index bowtie/${f}.sort.bam
	bamPEFragmentSize --maxFragmentLength 300 --binSize 1000 -p 4 -b bowtie/${f}.sort.bam -o FragmentSize/${f}_FragmentSize.svg
	samtools flagstat bowtie/${f}.sort.bam > log/${f}.txt
done

#Sample correlation checkup (spearman) by heatmap
multiBamSummary bins -bs 10000 -p 4 --verbose --bamfiles bowtie/H2AubEmb_Rep1.sort.bam bowtie/H2AubEmb_Rep2.sort.bam bowtie/H2Aub2Cells_Rep1.sort.bam bowtie/H2Aub2Cells_Rep2.sort.bam bowtie/H2AubSp_Rep1.sort.bam bowtie/H2AubSp_Rep2.sort.bam bowtie/H2AubEgg_Rep1.sort.bam bowtie/H2AubEgg_Rep2.sort.bam -o Correlation/ChipCorrelation.npz
plotCorrelation --removeOutliers --corData Correlation/ChipCorrelation.npz -c spearman -p heatmap --colorMap bwr --plotNumbers --labels H2AubEmb1 H2AubEmb2  H2Aub2Cells1 H2Aub2Cells2 H2AubSp1 H2AubSp2 H2AubOo1 H2AubOo2 -o Correlation/HeatmapCorChIP_Inputspearman.svg

# mark duplicates with Picard && index result bam
for f in `ls -1 bowtie/*.sort.bam | cut -d "/" -f 2 | sed 's/.sort.bam//'`
	do java -jar picard.jar MarkDuplicates I=bowtie/${f}.sort.bam O=bowtie/${f}.duplicates.bam M=bowtie/${f}.dup_metrics.txt
	   samtools index bowtie/${f}.duplicates.bam
done

# Remove PCR duplicates, then sort and index BAM
for f in `ls -1 bowtie/*.duplicates.bam  | cut -d "/" -f 2 | sed 's/.duplicates.bam//'`
	do samtools view -hb -F 0x400 -q 20 bowtie/${f}.duplicates.bam > bowtie/${f}.filter.bam
	   samtools sort -@ 6 bowtie/${f}.filter.bam -o bowtie/${f}.sort.filter.bam
	   samtools index bowtie/${f}.sort.filter.bam
	   samtools flagstat bowtie/${f}.sort.filter.bam > log/${f}.txt
done

# extracts reads from SNP calling : Genome 1 is maternal, Genome 2 is paternal 
for f in `ls -1 bowtie/*.sort.filter.bam  | cut -d "/" -f 2 | sed 's/.sort.filter.bam//'`
	do 
	SNPsplit --paired --snp_file SNPSplit/all_SNPs_PWK_PhJ_GRCm39.txt.gz -o PWK bowtie/${f}.sort.filter.bam
done


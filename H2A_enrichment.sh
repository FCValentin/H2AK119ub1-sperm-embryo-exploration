#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -q max-1m.q
#$ -e  ./log/
#$ -o ./log/

#3 is chr name, 4 begin of pos, 4+50 begin add length read, 1 is the name of read and abs($9) for fragment Length (TLEN) (distance 5' mate 1 3' mate 2 en paired end)
samtools view Merged/Input_Sp_WithoutDupl.sort.bam | awk '{print $3"\t"$4"\t"$4+50"\t"$1"\t"sqrt(($9)^2)}' > Merged/Input_Sp_WithoutDupl.sort.bed
samtools sort -n Merged/Input_Sp_WithoutDupl.sort.bam | samtools fixmate | bedtools bamtobed -bedpe > Merged/Input_Sp_WithoutDupl-bis.sort.bed
awk '{ if($5>135) print $0"\t"1; else if($5>30 && $5<=135) print $0"\t"2; else if($5<=30) print $0"\t"0}' Merged/Input_Sp_WithoutDupl.sort.bed > Merged/Input_Sp_WithoutDupl_flagNucl.bed
awk '{ if($5>95) print $0"\t"1; else if($5>30 && $5<=95) print $0"\t"2; else if($5<=30) print $0"\t"0}' Merged/Input_Sp_WithoutDupl.sort.bed > Merged/Input_Sp_WithoutDupl_H2A.bed

#selection of fragment length based on Fragment size output figure (BamPEFragmentsize)
awk '{ if($5<=90 && $5>50) print $0"\t"70; else if($5>140 && $5<165) print $0"\t"150; else if($5>100 && $5<=130) print $0"\t"110; else print $0"\t"0}' Merged/Input_Sp_WithoutDupl.sort.bed > Merged/Input_Sp_WithoutDupl_flag.bed

#split the flagged BEDPE file in 2 independant bed files depending H2A
nohup grep $'\t'1$ Merged/Input_Sp_WithoutDupl_flagNucl.bed | sort -k 1,1 -k2,2n > Merged/Input_Sp_WithoutDupl_Nucl.bed
nohup grep $'\t'2$ Merged/Input_Sp_WithoutDupl_flagNucl.bed | sort  -k 1,1 -k2,2n > Merged/Input_Sp_WithoutDupl_SubNucl.bed
nohup grep $'\t'1$ Merged/Input_Sp_WithoutDupl_H2A.bed | sort -k 1,1 -k2,2n > Merged/Input_Sp_WithoutDupl_H2APos.bed
nohup grep $'\t'2$ Merged/Input_Sp_WithoutDupl_H2A.bed | sort  -k 1,1 -k2,2n > Merged/Input_Sp_WithoutDupl_H2ANeg.bed

#coverage of fragment length reads on each genome bins
bedtools intersect -sorted -a Genome/XLaevis_250bpbin_slid50.bed -b Merged/Input_Sp_WithoutDupl_Nucl.bed -c > Fragment/Input_Sp_WithoutDupl_Nucl_slid50.bed
bedtools intersect -sorted -a Genome/XLaevis_250bpbin_slid50.bed -b Merged/Input_Sp_WithoutDupl_SubNucl.bed -c > Fragment/Input_Sp_WithoutDupl_SubNucl_slid50.bed
bedtools intersect -sorted -a Genome/XLaevis_250bpbin_slid50.bed -b Merged/Input_Sp_WithoutDupl_H2APos.bed -c > Fragment/Input_Sp_WithoutDupl_H2APos_slid50.bed
bedtools intersect -sorted -a Genome/XLaevis_250bpbin_slid50.bed -b Merged/Input_Sp_WithoutDupl_H2ANeg.bed -c > Fragment/Input_Sp_WithoutDupl_H2ANeg_slid50.bed

#Format setup
paste Fragment/Input_Sp_WithoutDupl_Nucl_slid50.bed Fragment/Input_Sp_WithoutDupl_SubNucl_slid50.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"($4+$8)}' > Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_Nucl_vs_SemiNucl.full.bedgraph
awk '{if($6>=5) print $0}' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_Nucl_vs_SemiNucl.full.bedgraph > Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.full.bedgraph
paste Fragment/Input_Sp_WithoutDupl_H2APos_slid50.bed Fragment/Input_Sp_WithoutDupl_H2ANeg_slid50.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"($4+$8)}' > Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_H2A_vs_NonH2A.full.bedgraph
awk '{if($6>=5) print $0}' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_H2A_vs_NoH2A.full.bedgraph > Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.full.bedgraph

############################################################################################
##### R Script PValues : Output a .probabilities file : check ProcessingChIPdata.R #########
############################################################################################

# P-Value discretisation
awk '
    {
    if ($6>0) {p1=0;p2=0}
    if($6>0 && $9>0.1) {p1=0;}
    if($6>0 && $9<=0.1 && $9>0.05) {p1=1;} 
    if($6>0 && $9<=0.05 && $9>0.001) {p1=2;}
    if($6>0 && $9<=0.001){p1=3;}
    if($6>0 && $10>0.1) {p2=0;}
    if($6>0 && $10<=0.1 && $10>0.05) {p2=1;} 
    if($6>0 && $10<=0.05 && $10>0.001) {p2=2;}
    if($6>0 && $10<=0.001){p2=3;}
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"p1"\t"p2}
' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities >  Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities.discrete

awk '
    {
    if ($6>0) {p1=0;p2=0}
    if($6>0 && $9>0.1) {p1=0;}
    if($6>0 && $9<=0.1 && $9>0.05) {p1=1;} 
    if($6>0 && $9<=0.05 && $9>0.001) {p1=2;}
    if($6>0 && $9<=0.001){p1=3;}
	if($6>0 && $9<=0.0001 && $9>0.00001){p1=4;}
	if($6>0 && $9<=0.00001){p1=5;}
    if($6>0 && $10>0.1) {p2=0;}
    if($6>0 && $10<=0.1 && $10>0.05) {p2=1;} 
    if($6>0 && $10<=0.05 && $10>0.001) {p2=2;}
    if($6>0 && $10<=0.001 && $10>0.0001){p2=3;}
	if($6>0 && $10<=0.0001 && $10>0.00001){p2=4;}
	if($6>0 && $10<=0.00001){p2=5;}
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"p1"\t"p2}
' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities >  Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities.discrete

###############################################################
##### R Script HMM Model : check ProcessingChIPdata.R #########
###############################################################

#Merged bins that overlap between them into a unique one, at different discrete p-value level 
awk '{if($11>0) print $0}' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities.discrete | mergeBed > IGV/H2Aenriched_0.1.bed
awk '{if($12>0) print $0}' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities.discrete | mergeBed > IGV/H2Adepleted_0.1.bed
awk '{if($11>0) print $0}' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities.discrete | mergeBed > IGV/NucleosomeEnriched_0.1.bed
awk '{if($12>0) print $0}' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities.discrete | mergeBed > IGV/NucleosomeDepleted_0.1.bed

awk '{if($11>1) print $0}' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities.discrete | mergeBed > IGV/H2Aenriched_0.05.bed
awk '{if($12>1) print $0}' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities.discrete | mergeBed > IGV/H2Adepleted_0.05.bed
awk '{if($11>1) print $0}' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities.discrete | mergeBed > IGV/NucleosomeEnriched_0.05.bed
awk '{if($12>1) print $0}' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities.discrete | mergeBed > IGV/NucleosomeDepleted_0.05.bed

awk '{if($11>2) print $0}' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities.discrete | mergeBed > IGV/H2Aenriched_0.001.bed
awk '{if($12>2) print $0}' Fragment/Input.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NoH2A.probabilities.discrete | mergeBed > IGV/H2Adepleted_0.001.bed
awk '{if($11>2) print $0}' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities.discrete | mergeBed > IGV/NucleosomeEnriched_0.001.bed
awk '{if($12>2) print $0}' Fragment/Input.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.probabilities.discrete | mergeBed > IGV/NucleosomeDepleted_0.001.bed

#Compute genome coverage of each particle enrichment
for f in `ls -1 IGV/Nucl*_*.bed | cut -d "/" -f 2 | sed 's/.bed//'`
	do bedtools coverage -a ../GFF/XL9_Chromosomes.bed -b IGV/${f}.bed > Coverage/Coverage_${f}.bed
done

#Nucleosome or H2A enrichment on H2AK119ub1 peaks
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b NucleosomeEnriched_0.001.bed > NucleosomeEnriched_0.001.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b NucleosomeEnriched_0.05.bed > NucleosomeEnriched_0.05.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b NucleosomeEnriched_0.1.bed > NucleosomeEnriched_0.1.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b NucleosomeDepleted_0.001.bed > NucleosomeDepleted_0.001.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b NucleosomeDepleted_0.05.bed > NucleosomeDepleted_0.05.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b NucleosomeDepleted_0.1.bed > NucleosomeDepleted_0.1.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b H2Aenriched_0.001.bed > H2Aenriched_0.001.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b H2Aenriched_0.05.bed > H2Aenriched_0.05.coverage
bedtools coverage -a Merged_Untreated_Sp_H2Aub_ChIP_001.XL92_peaks.bed -b H2Aenriched_0.1.bed > H2Aenriched_0.1.coverage

#sperm H2AK119ub1 signal on Nucleosome or H2A enrichment 
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeEnriched_0.001.bed > FragmentSize/Input_NucleosomeEnriched_0.001.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeEnriched_0.05.bed > FragmentSize/Input_NucleosomeEnriched_0.05.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeEnriched_0.1.bed > FragmentSize/Input_NucleosomeEnriched_0.1.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeDepleted_0.001.bed > FragmentSize/Input_NucleosomeDepleted_0.001.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeDepleted_0.05.bed > FragmentSize/Input_NucleosomeDepleted_0.05.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeDepleted_0.1.bed > FragmentSize/Input_NucleosomeDepleted_0.1.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/H2Aenriched_0.001.bed > FragmentSize/Input_H2Aenriched_0.001.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/H2Aenriched_0.05.bed > FragmentSize/Input_H2Aenriched_0.05.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b IGV/H2Aenriched_0.1.bed > FragmentSize/Input_H2Aenriched_0.1.bam
bedtools intersect -f 0.5 -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeEnriched_0.001.bed > FragmentSize/Input_NucleosomeEnriched_0.001_MidOverlap.bam
bedtools intersect -f 0.5 -a Merged/Input_Sp_WithoutDupl.bam -b IGV/NucleosomeDepleted_0.001.bed > FragmentSize/Input_NucleosomeDepleted_0.001_MidOverlap.bam
bedtools intersect -f 0.5 -a Merged/Input_Sp_WithoutDupl.bam -b IGV/H2Aenriched_0.001.bed > FragmentSize/Input_H2Aenriched_0.001_MidOverlap.bam

#Check fragment size of isolated reads on fragment size (control test)
for f in `ls -1 FragmentSize/*_MidOverlap.bam  | cut -d "/" -f 2 | sed 's/.bam//'`
	do 
	   bamCoverage -bs 50 -v -p 4 -b FragmentSize/${f}.bam -o IGV/${f}.bw
	   samtools sort -@ 6 FragmentSize/${f}.bam -o FragmentSize/${f}.sort.bam
	   samtools index FragmentSize/${f}.sort.bam
	   bamPEFragmentSize -p 8 --plotFileFormat pdf -b FragmentSize/${f}.sort.bam -T "Fragment size of PE RNA-seq datas" --maxFragmentLength 300 -o FragmentSize/${f}_Hist.pdf --table FragmentSize/${f}.tsv --samplesLabel ${f}
done

#extract Sperm reads fragment size reads on H2A enriched or subNucleosome genomic enriched regions
sort -k1,1 -k2,2n Merged/Input_Sp_WithoutDupl_SubNucl.bed | mergeBed > Merged/Input_Sp_WithoutDupl_SubNucl.sorted.bed 
sort -k1,1 -k2,2n Merged/Input_Sp_WithoutDupl_H2APos.bed | mergeBed > Merged/Input_Sp_WithoutDupl_H2APos.sorted.bed 
sort -k1,1 -k2,2n Merged/Input_Sp_WithoutDupl_H2ANeg.bed | mergeBed > Merged/Input_Sp_WithoutDupl_H2ANeg.sorted.bed 
bedtools intersect -a Merged/Input_Sp_WithoutDupl.bam -b Merged/Input_Sp_WithoutDupl_Nucl.bed > IGV/Input_Nucleosome.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.sort.bam -b Merged/Input_Sp_WithoutDupl_SubNucl.sorted.bed > IGV/Input_SubNucleosome.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.sort.bam -b Merged/Input_Sp_WithoutDupl_H2APos.sorted.bed > IGV/Input_WithH2A.bam
bedtools intersect -a Merged/Input_Sp_WithoutDupl.sort.bam -b Merged/Input_Sp_WithoutDupl_H2ANeg.sorted.bed > IGV/Input_WithoutH2A.bam
samtools sort -@ 6 IGV/Input_Nucleosome.bam -o IGV/Input_Nucleosome.sort.bam
samtools sort -@ 6 IGV/Input_SubNucleosome.bam -o IGV/Input_SubNucleosome.sort.bam
samtools sort -@ 6 IGV/Input_WithH2A.bam -o IGV/Input_WithH2A.sort.bam
samtools sort -@ 6 IGV/Input_WithoutH2A.bam -o IGV/Input_WithoutH2A.sort.bam

#Compute H2A and Nucleosome enrichment IGV track 
samtools index IGV/Input_Nucleosome.sort.bam
samtools index IGV/Input_SubNucleosome.sort.bam
bamCompare --binSize 50 -p 4 -v -of bigwig -b1 IGV/Input_Nucleosome.sort.bam -b2 IGV/Input_SubNucleosome.sort.bam -o IGV/NuclPos_vs_SubNuclNeg_log2ratio.bw

samtools index IGV/Input_WithH2A.sort.bam
samtools index IGV/Input_WithoutH2A.sort.bam
bamCompare --binSize 50 -p 4 -v -of bigwig -b1 IGV/Input_WithH2A.sort.bam -b2 IGV/Input_WithoutH2A.sort.bam -o IGV/H2APos_vs_H2ANeg_log2ratio.bw


###########################################
##### DNA Ladder ratio ChIP/Input #########
###########################################

#RATIO H2Aub => 10807/484 = 22.33 (from 601_DNA reads data extracted, ChIP/Input 601_DNA reads)

#extract ChIP and Input sperm reads on Xenopus laevis bins
awk '{print $1"\t"($2+100)"\t"($3-100)}' DNALadder/XLaevis_250bpbin_slid50.bed > DNALadder/XLaevis_50bpbin.bed

#For each sample file, compute Ladder normalized signal from ChIP. Be careful to edit RATIO H2Aub correction, at 48400/10807 in that case
for f in `ls -1 DNALadder/*_INPUT_Input.XL92.sort.filter.bam | cut -d "/" -f 2 | sed 's/_INPUT_Input.XL92.sort.filter.bam//'`
	do
	bedtools intersect -a DNALadder/XLaevis_250bpbin.bed -b DNALadder/${f}_INPUT_Input.XL92.sort.filter.bam -c > DNALadder/${f}_INPUT_Input.GenomeCov.bed
	bedtools intersect -a DNALadder/XLaevis_250bpbin.bed -b DNALadder/${f}_H2Aub_ChIP.XL92.sort.filter.bam -c > DNALadder/${f}_H2Aub_ChIP.GenomeCov.bed
	paste DNALadder/${f}_H2Aub_ChIP.GenomeCov.bed DNALadder/${f}_INPUT_Input.GenomeCov.bed | awk '{if($8>0){print $1"\t"$2"\t"$3"\t"$4"\t"($8)*(2)"\t"($4/(($8)*(2)))"\t"($4/(($8)*(2)))*(48400/10807);}if($8=0){print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t""0"}}' > DNALadder/${f}_LadderCov.bed
	awk '{print $1"\t"$2"\t"$3"\t"$7}' DNALadder/${f}_LadderCov.bed > DNALadder/${f}_LadderCov.bg
	./bedGraphtoBigWig DNALadder/${f}_LadderCov.bg chromSizesXL9.bed DNALadder/${f}_LadderCov.bw
done

#Ladder
computeMatrix reference-point -p 4 --referencePoint center --verbose -S cutadapt_oo-WT_LadderCov.bw cutadapt_Unt_LadderCov.bw cutadapt_oo-U21_LadderCov.bw cutadapt_Untreated_Sp_LadderCov.bw cutadapt_Egg-oo-WT_LadderCov.bw cutadapt_Egg-Unt_LadderCov.bw cutadapt_Egg-oo-U21_LadderCov.bw cutadapt_EggExtract_Sp_LadderCov.bw -R USP21sensitivesTSS.bed MZTSS.bed MaternalOnlyTSS.bed ZygoticTSS.bed OthersTSS.bed GenesTSS.bed --binSize 50 --missingDataAsZero --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --skipZeros -out Droso_USP21Sensitive_LadderH2Aub.mat.gz --outFileNameMatrix Droso_USP21Sensitive_LadderH2Aub.tab --outFileSortedRegions Droso_USP21Sensitive_LadderH2Aub.region.bed
plotProfile --samplesLabel oo-WTH2Aub oo-UntH2Aub oo-USP21H2Aub oldoo-UntH2Aub egg-WTH2Aub egg-UntH2Aub egg-USP21H2Aub oldegg-UntH2Aub --regionsLabel USP21_DEGenes MZGenes MaternalNonZygoticGenes ZygoticGenes OthersGenes AllGenes -m Droso_USP21Sensitive_LadderH2Aub.mat.gz -out Droso_USP21Sensitive_LadderH2Aub.profile.svg --outFileNameData Droso_USP21Sensitive_LadderH2Aub.tsv

#H2AK119ub1 calibration
paste DNALadder/cutadapt_Untreated_Sp_H2Aub_ChIP.250bp.bed DNALadder/cutadapt_Untreated_Sp_INPUT_Input.250bp.bed | awk '{if($8>0){print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"($4/($8))"\t"($4/($8))*(48400/10807);}if($8=0){print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t""0"}}' > DNALadder.bed
awk '{print $1"\t"$2"\t"$3"\t"$7}' DNALadder.bed  > DNALadder.bg
awk '{if($4>100){print $1"\t"$2"\t"$3"\t""100"}if($4<=100){print $1"\t"$2"\t"$3"\t"$4}}' DNALadder.bg  > DNALadder_MAX.bg
./bedGraphtoBigWig DNALadder_MAX.bg chromSizesXL9.bed DNALadder.bw

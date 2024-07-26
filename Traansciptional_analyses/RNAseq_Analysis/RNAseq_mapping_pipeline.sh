## RNA-seq Pipeline to Map with BT2 Df, BT2 k100, and BT2 k100 (51mer Filetered) for 0-12h embryos ###
### Same pipeline used to map RNA-seq brains data ###
######################################################################################################

module load cutadapt
module load fastqc
module load bowtie2/2.3.5.1
module load samtools/1.16.1
module load bedtools
module load ucsc_genome 

REF=dmel_scaffold2_plus0310.fasta
################## PRE-PROCESSING ##################
fastqc ${INpath}${INfileprefix}.fastq.gz
cutadapt -j 24 -q 20 -m 100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${INpath}${INfileprefix}_R1_trimmed.fastq -p $ {INpath}${INfileprefix}_R2_trimmed.fastq ${INpath}${INfileprefix}_R1.fastq ${INpath}${INfileprefix}_R2.fastq
fastqc ${INpath}${INfileprefix}_trimmed.fastq

#### Bowtie2 Df Mapping ####
############################
bowtie2 -p 24 -t -x $REF -1 Emb_RNAseq_R1_trimmed.fastq.gz -2 Emb_RNAseq_R2_trimmed.fastq.gz | \
        samtools view -@ 24 -F 1548 -bS | samtools sort -O bam -o \
        RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_sort.bam

#### Bowtie2 K-100 Mapping ####
###############################
bowtie2 -p 24 -t -k 100 -x $REF -1 Emb_RNAseq_R1_trimmed.fastq.gz -2 Emb_RNAseq_R2_trimmed.fastq.gz | \
        samtools view -@ 24 -F 1548 -bS | samtools sort -O bam -o \
        RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_sort.bam

#### Split Bt2 Df Files ####
############################
bedtools genomecov -bg -split -strand + -ibam RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_sort.bam contigs.size RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_plus.bedgraph
bedtools genomecov -bg -split -strand - -ibam RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_sort.bam contigs.size RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_m.bedgraph
awk '{ $4=$4*-1; print }' RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_m.bedgraph > RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_minus.bedgraph

export LC_COLLATE=C
sort -k1,1 -k2,2n RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_plus.bedgraph > RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_plus_sort.bedgraph
sort -k1,1 -k2,2n RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_minus.bedgraph > RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_minus_sort.bedgraph
bedGraphToBigWig RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_plus_sort.bedgraph contigs.size RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_plus.bigwig
bedGraphToBigWig RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_minus_sort.bedgraph contigs.size RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_minus.bigwig

#### Split Bt2 K-100 Files ####
###############################
bedtools genomecov -bg -split -strand + -ibam RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_sort.bam contigs.size RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_plus.bedgraph
bedtools genomecov -bg -split -strand - -ibam RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_sort.bam contigs.size RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_m.bedgraph
awk '{ $4=$4*-1; print }' RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_m.bedgraph > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_minus.bedgraph

export LC_COLLATE=C
sort -k1,1 -k2,2n RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_plus.bedgraph > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_plus_sort.bedgraph
sort -k1,1 -k2,2n RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_minus.bedgraph > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_minus_sort.bedgraph
bedGraphToBigWig RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_plus_sort.bedgraph contigs.size RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_plus.bigwig
bedGraphToBigWig RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_minus_sort.bedgraph contigs.size RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_minus.bigwig

### Convert bam to bed file ###
bedtools bamtobed -i RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548_sort.bam > RNAseq_Emb_q20-m100_BT2-Df-Dmel-F1548.bed
bedtools bamtobed -i RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_sort.bam > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548.bed

#### Overlap 51-mer Meryl ####
##############################
overlapSelect -overlapBases=51 dmel_assembly_51mer_SINGLEv1.3.meryl_merge.bed RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548.bed \
        RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge.bed

#### Split K100 filtered Files ####
###################################
awk '$6 == "+"' RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge.bed | genomeCoverageBed -i contigs.size > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_plus.bedgraph
awk '$6 == "-"' RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge.bed | genomeCoverageBed -i contigs.size > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_m.bedgraph
awk '{ $4=$4*-1; print }' RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_m.bedgraph > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_minus.bedgraph

export LC_COLLATE=C
sort -k1,1 -k2,2n RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_plus.bedgraph > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_plus_sort.bedgraph
sort -k1,1 -k2,2n RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_minus.bedgraph > RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_minus_sort.bedgraph

bedGraphToBigWig RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_plus_sort.bedgraph contigs.size \
        RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_plus.bigwig
bedGraphToBigWig RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_minus_sort.bedgraph contigs.size \
        RNAseq_Emb_q20-m100_BT2-K100-Dmel-F1548_OVER_Meryl_uniq51mer_Dmel-merge_minus.bigwig

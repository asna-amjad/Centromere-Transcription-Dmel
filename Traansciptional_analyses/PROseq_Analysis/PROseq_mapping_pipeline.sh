## RNA-seq Pipeline to Map with BT2 Df, BT2 k100, and BT2 k100 (51mer Filetered) for 0-12h embryos ###
### Same pipeline used to map RNA-seq brains data ###
######################################################################################################

module load cutadapt
module load fastqc
module load bowtie2/2.3.5.1
module load samtools/1.16.1
module load bedtools
module load ucsc_genome

mkdir fastq_Trim

INpath=/path/to/input/directory/
INfileprefix=SAMPLE1
DmelanogasterIndex=/labs/Mellone/Asna/Assembly/bowtie2_Index/dmel_scaffold2_plus0310.fasta
OUTpath=/path/to/output/directory/

################## Bowtie2 Default Mapping Pipeline ##################
######################################################################

### PRE-PROCESSING ###
fastqc ${INpath}${INfileprefix}.fastq
cutadapt -j 24 -q 20 -m 20 -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC -o ${INpath}${INfileprefix}_R1_trim.fastq.gz -p $ {INpath}${INfileprefix}_R2_Trim.fastq.gz \
  ${INpath}${INfileprefix}_R1.fastq.gz ${INpath}${INfileprefix}_R2.fastq.gz

### FILTERING OUT HeLa Cells SPIKE-INS IN PROSEQ DATA ONLY ###
bowtie2 -p 24 -t -x Dmel_hg38_combine -1 fastq_Trim/ProSeq_5_R1_Trim.fastq.gz \
       -2 fastq_Trim/ProSeq_5_R2_Trim.fastq.gz | \
       samtools view -@ 24 -F 1548 -bS | \
       samtools sort -@ 24 -m 3G -O bam -o PROseq_Emb_BT2_Df_CombineHgDmel_sort.bam

(samtools view -@ 24 -H PROseq_Emb_BT2_Df_CombineHgDmel_sort.bam
samtools view -@ 24 PROseq_Emb_BT2_Df_CombineHgDmel_sort.bam | grep -f contigs.txt) | \
samtools view -@ 24 -b > PROseq_Emb_BT2_Df_Dmel_sort.bam

## Step 1: Get concordantly aligned reads
samtools view -@ 24 -hf 0x2 PROseq_Emb_BT2_Df_Dmel_sort.bam |
awk '($1 ~/^@/ || $2==99 || $2==147)' |
samtools sort -n -@ 24 |
bedtools bamtobed -bedpe -i stdin |
awk '{OFS="\t"} {print $1,$2,$6,$7,$8,"-"}' > PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk.bed
samtools view -@ 24 -hf 0x2 PROseq_Emb_BT2_Df_Dmel_sort.bam |
awk '($1 ~/^@/ || $2==83 || $2==163)' |
samtools sort -n -@ 24 |
bedtools bamtobed -bedpe -i stdin |
awk '{OFS="\t"} {print $1,$2,$6,$7,$8,"+"}' >> PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk.bed

Step 2: Get 3' ends per pair (position of polymerase) & convert to bedgraph
cat PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk.bed |
    awk '{OFS="\t"} $6=="+" {print $1,($3-1),$3,$4,$5"_"($3-$2),$6}; $6=="-" {print $1,$2,($2+1),$4,$5"_"($3-$2),$6}' |
    sort --parallel 12 -k1,1 -k2,2n > PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends.bed
##Convert to bedgraph:
awk '$6 == "+"' PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs.size > \
    PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_plus.bedgraph
awk '$6 == "-"' PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs.size > \
    PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_m.bedgraph
awk '{OFS="\t"}{ $4=$4*-1; print }' PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_m.bedgraph > \
    PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_minus.bedgraph

export LC_COLLATE=C
    sort -k1,1 -k2,2n PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_plus.bedgraph > PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_plus_sort.bedgraph
    sort -k1,1 -k2,2n PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_minus.bedgraph > PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_minus_sort.bedgraph
bedGraphToBigWig PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_plus_sort.bedgraph contigs.size PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_plus.bigwig
bedGraphToBigWig PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_minus_sort.bedgraph contigs.size PROseq_Emb_BT2_Df_Dmel_nameSort_withAwk_3ends_minus.bigwig  

### Bowtie2 K100 Pipeline ###
#############################
### Reverse Complement
seqkit seq --reverse --complement ProSeq_5_R1_Trim.fastq.gz > PROseq_Trim_R1_RC.fastq.gz
### Align to Hg28 genome 
bowtie2 -p 24 -t --very-sensitive -x fastq_Trim/AlignProSeqDec23/hg38_bowtie2/Human_genome \
       -U PROseq_Trim_R1_RC.fastq.gz | samtools view -@ 24 -bS | \
       samtools sort -@ 24 -m 3G -O bam -o PROseq_Trim_R1_RC_BT2-vs-Hg_sort.bam

### Pull out UNMAPPED PROseq reads only from sort.bam and convert to fastq
samtools view -@ 24 -b -f 4 -o PROseq_Trim_R1_RC_BT2-vs-Hg_sort_unmapped.bam PROseq_Trim_R1_RC_BT2-vs-Hg_sort.bam
samtools fastq -@ 24 PROseq_Trim_R1_RC_BT2-vs-Hg_sort_unmapped.bam > PROseq_Trim_R1_RC_BT2-vs-Hg_sort_unmapped.fastq

### Map to Drosophila
bowtie --threads 24 -k 100 --sam -q -x Bt1_Index/dmel_het_assembly.fasta \
       PROseq_Trim_R1_RC_BT2-vs-Hg_sort_unmapped.fastq | samtools view -@ 24 -bS | samtools sort -m 3G -O bam -o PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort.bam

bedtools bamtobed -i PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort.bam > PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM.bed
sort -k1,1 -k2,2n PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM.bed > PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort.bed

create 3'-end BEDs (FOR KMER OVERLAPS & getting counts - add to directory)
awk '{OFS="\t"} ($6=="+") {print $1,($3-1),$3,".",($3-$2),$6}; ($6=="-") {print $1,$2,($2+1),".",($3-$2),$6}' PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort.bed |
        sort -k1,1 -k2,2n > PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime.bed

#generate bedgraphs (for tracks)
awk '$6 == "+"' PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs.size > \
       PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_plus.bedgraph
awk '$6 == "-"' PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs.size > \
       PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_m.bedgraph
awk '{OFS="\t"}{ $4=$4*-1; print }' PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_m.bedgraph > PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_minus.bedgraph

export LC_COLLATE=C
sort -k1,1 -k2,2n PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_plus.bedgraph > PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_plus_sort.bedgraph
sort -k1,1 -k2,2n PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_minus.bedgraph > PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_minus_sort.bedgraph

bedGraphToBigWig PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_plus_sort.bedgraph contigs.size PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_plus.bigwig
bedGraphToBigWig PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_minus_sort.bedgraph contigs.size PROseq_Trim_R1_RC_HumanFilt_BT-k100-DM_sort_3prime_minus.bigwig


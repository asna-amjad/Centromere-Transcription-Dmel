#!/bin/bash
#SBATCH --job-name=strand
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=asna.amjad@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

module load samtools/1.9 
module load bedtools
module load ucsc_genome
module load GenomeBrowser/20180626 

mkdir -p ~/.tmp
export TMPDIR=~/.tmp

## 1) Start with concordantly aligned minus strand read pairs
## ============================================================

samtools view -hf 0x2 ../jockey-3_hybcap2dmel_q30.bam |
  awk '($1 ~/^@/ || $2==99 || $2==147)' |
  samtools sort -n -@ 12 |
  bedtools bamtobed -bedpe -i stdin |
  awk '{OFS="\t"} {print $1,$2,$6,$7,$8,"-"}' > jockey-3_hybcap2dmel_q30.bed

## 2) Add in concordantly aligned plus strand read pairs
## ======================================================

samtools view -hf 0x2 ProSeq_embryo_Iso1_bowtie2_Df_sorted.bam |
  awk '($1 ~/^@/ || $2==83 || $2==163)' |
  samtools sort -n -@ 12 |
  bedtools bamtobed -bedpe -i stdin |
  awk '{OFS="\t"} {print $1,$2,$6,$7,$8,"+"}' >> ProSeq_embryo_Iso1_bowtie2_Df_nameSort.bed


## 3) 3' ends of alignments
## ==========================

cat ../jockey-3_hybcap2dmel.bed |
  awk '{OFS="\t"} $6=="+" {print $1,($3-1),$3,$4,$5"_"($3-$2),$6}; $6=="-" {print $1,$2,($2+1),$4,$5"_"($3-$2),$6}' |
  sort -k1,1 -k2,2n > jockey-3_hybcap2dmel_sort.bed

## 4) make bedgraph of 3' ends (position of polymerase)
## =====================================================

## ***sort genome file first to avoid chaos downstream***
## -------------------------------------------------------

#export LC_COLLATE=C
#sort -k1,1 -k2,2n contigs.size > contigs_sort.size

## *** plus strand ***
## --------------------
awk '$6 == "+"' ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs_sort.size > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bedgraph

## *** minus strand ***
## --------------------
awk '$6 == "-"' ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs_sort.size > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m.bedgraph

## *** negate minus strand ***
## ---------------------------
awk '{ $4=$4*-1; print }' ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bedgraph

## 5) sort 3' end bedgraphs
## =========================

#export LC_COLLATE=C
sort  -k1,1 -k2,2n ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus_sort.bedgraph
sort  -k1,1 -k2,2n ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus_sort.bedgraph

## (non-negated minus strand) 
sort -k1,1 -k2,2n ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m_sort.bedgraph

## 6) bedgraph to bigwig
## ======================


bedGraphToBigWig ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus_sort.bedgraph contigs_sort.size ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig
bedGraphToBigWig ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus_sort.bedgraph contigs_sort.size ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig


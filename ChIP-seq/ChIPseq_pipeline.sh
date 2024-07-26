#!/bin/bash
#SBATCH --job-name=newAlign
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=150G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=asna.amjad@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#module load TrimGalore
#module load cutadapt
#module load fastqc
module load bowtie2/2.3.5.1
module load samtools/1.16.1
module load bedtools
#module load ucsc_genome
#module load GenomeBrowser/20180626 

mkdir -p ~/.tmp
export TMPDIR=~/.tmp

mkdir fastq_Trim

trim_galore -gzip --paired --illumina -length 35 --phred33 -o fastq_Trim SRR7588744/SRR7588744_1.fastq SRR7588744/SRR7588744_2.fastq
trim_galore -gzip --paired --illumina -length 35 --phred33 -o fastq_Trim SRR7588744/SRR7588745_1.fastq SRR7588744/SRR7588745_2.fastq
 
mkdir Mapping
bowtie2 -p 48 -t -x /labs/Mellone/Asna/Assembly/bt2_Index/dmel_scaffold2_plus0310.fasta \
        -1 SRR7588745_1_val_1.fq.gz  \
        -2 SRR7588745_2_val_2.fq.gz | samtools view -@ 48 -bS -f 0x2 | \
        samtools sort -@ 48 -n -m 3G -O bam -o Mapping/CenpA_ChipSeq_BT2-vs-Dmel_nameSort.bam

bowtie2 -p 48 -t -x /labs/Mellone/Asna/Assembly/bt2_Index/dmel_scaffold2_plus0310.fasta \
       -1 /home/FCAM/aamjad/fly_T2T_assemblies/CENP-A_OreR_ChIP-seq/fastq_Trim/SRR7588744_1_val_1.fq.gz  \
       -2 /home/FCAM/aamjad/fly_T2T_assemblies/CENP-A_OreR_ChIP-seq/fastq_Trim/SRR7588744_2_val_2.fq.gz | samtools view -@ 24 -bS -f 0x2 | \
       samtools sort -@ 48 -n -m 3G -O bam -o Mapping/Input_ChipSeq_BT2-vs-Dmel_nameSort.bam

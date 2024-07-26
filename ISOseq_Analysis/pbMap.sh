#!/bin/bash
#SBATCH --job-name=map
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=asna.amjad@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

module load pb-assembly

pbmm2 align -j 8 --preset ISOSEQ --sort --unmapped dmel_scaffold2_plus0310.mmi /labs/Mellone/Asna/ISOseq/FASTQ_data_files/jockey-3_hybcap_flnc.fa jockey-3_hybcap_flnc.aligned.bam --log-level INFO

# Centromere-Transcription-Dmel

This repository is for all the analyses related to my project of exploring centromere transcription in Drosophila melanogaster datasets (PROseq, RNAseq, and CUT&Tag). All analyses are mainly done on the University of Connecticut HPC cluster under SLURM software.

Each code block is uploaded as a basic set of commands, and was run on the University of Connecticut HPC cluster under the SLURM management software. This code block is included here: 
#!/bin/bash 
#SBATCH --job-name=X 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 15 
#SBATCH --partition=general 
#SBATCH --qos=general 
#SBATCH --mem=100g 
#SBATCH -o %x_%j.out 
#SBATCH -e %x_%j.err 

This block was modified for memory and performance requirements accordingly.

Meryl-associated code and unique kmer: https://github.com/bmellone/Dmel-Centromere-Transcription/tree/main/Unique_Kmers

Drosophila heterochromatin assembly (Chang et al. PLoS Bio, 2019 DOI:https://doi.org/10.1371/journal.pbio.3000241) doi:10.5061/dryad.rb1bt3j/File.S8.Chang_et_al.fasta Drosophila heterochromatin repeat annotations (Chang et al. PLoS Bio, 2019 DOI:https://doi.org/10.1371/journal.pbio.3000241) doi:10.5061/dryad.rb1bt3j/File.S9.Chang_et_al.gff.txt

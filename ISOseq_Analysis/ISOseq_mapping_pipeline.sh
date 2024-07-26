#!/bin/bash
#SBATCH --job-name=bw
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

#module load minimap2
module load bedtools
module load deeptools
module load ucsc_genome

chromInfo=~/assembly/Het_chromatin_assembly/contigs.size 

##### Map to Assembly with MiniMap2 #####
#########################################
minimap2 -t 12 -ax map-hifi --MD /home/FCAM/aamjad/assembly/dmel_scaffold2_plus0310.fasta ../jockey-3_hybcap_flnc.fa -o jockey-3_hybcap_flnc_minimap2.sam

samtools view -@ 12 -bq 30 -F 2308 jockey-3_hybcap_flnc_minimap2.sam -o jockey-3_hybcap_flnc_minimap2.bam
samtools index -@ 12 jockey-3_hybcap_flnc_minimap2.bam

#minimap2 -t 12 -ax map-hifi --MD /labs/Mellone/Asna/ISOseq/FASTQ_data_files/J3_consensus_promoter.fasta /labs/Mellone/Asna/ISOseq/FASTQ_data_files/jockey-3_hybcap_flnc.fa \
#        -o jockey3_hybcap_promoter.sam

#bamCoverage -b jockey-3_hybcap_flnc_minimap2_sorted.bam -bs 1 --filterRNAstrand forward -o jockey-3_hybcap_flnc_minimap2_fwd.bw
#bamCoverage -b jockey-3_hybcap_flnc_minimap2_sorted.bam -bs 1 --filterRNAstrand reverse -o jockey-3_hybcap_flnc_minimap2_rev.bw

###############################################################
#Convert sort.bam to bedgraph, sort, then to bigwig for viewing:
    ##Versions: bedtools/2.29.0, GenomeBrowser/20180626
    ##bedgraphtobigwig documentation: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
    ##PROseq commands:

bedtools bamtobed -i jockey-3_hybcap_flnc_minimap2_sorted.bam > Coverage/jockey-3_hybcap_flnc_minimap2_sorted.bed

bedtools genomecov -bg -ibam jockey-3_hybcap_flnc_minimap2_sorted.bam > Coverage/jockey-3_hybcap_flnc_minimap2.bedgraph
export LC_COLLATE=C
sort -k1,1 -k2,2n Coverage/jockey-3_hybcap_flnc_minimap2.bedgraph > Coverage/jockey-3_hybcap_flnc_minimap2_sort.bedgraph
             
bedGraphToBigWig align/PROseq_RevComp_BT2-vs-dM_BT2-k100-F_1548_sort.bedgraph ~/assembly/contigs.size \
       align/bigWig/PROseq_RevComp_BT2-vs-dM_BT2-k100-F_1548_sort.bigwig

genomeCoverageBed -bg -strand + -i Coverage/jockey-3_hybcap_flnc_minimap2.bed -g ${chromInfo} > Coverage/jockey-3_hybcap_flnc_minimap2_plus.bdg

genomeCoverageBed -bg -strand - -i Coverage/jockey-3_hybcap_flnc_minimap2.bed -g ${chromInfo} > Coverage/jockey-3_hybcap_flnc_minimap2_minus.bdg


awk '{OFS="\t"}{$NF *= -1; print}' Coverage/jockey-3_hybcap_flnc_minimap2_minus.bdg > Coverage/jockey-3_hybcap_flnc_minimap2_minus_mod.bdg

export LC_COLLATE=C
sort -k1,1 -k2,2n Coverage/jockey-3_hybcap_flnc_minimap2_plus.bdg > Coverage/jockey-3_hybcap_flnc_minimap2_plus_sorted.bdg
sort -k1,1 -k2,2n Coverage/jockey-3_hybcap_flnc_minimap2_minus_mod.bdg > Coverage/jockey-3_hybcap_flnc_minimap2_minus_sorted.bdg

bedGraphToBigWig Coverage/jockey-3_hybcap_flnc_minimap2_plus_sorted.bdg ${chromInfo} Coverage/jockey-3_hybcap_flnc_minimap2_plus_sorted.bw
bedGraphToBigWig Coverage/jockey-3_hybcap_flnc_minimap2_minus_sorted.bdg ${chromInfo} Coverage/jockey-3_hybcap_flnc_minimap2_minus_sorted.bw 

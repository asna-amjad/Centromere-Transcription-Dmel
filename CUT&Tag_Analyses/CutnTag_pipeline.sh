## CUT&Tag Pipeline to Map with Bowtie2 for 0-12h embryos ###
#############################################################

module load TrimGalore/0.6.7 cutadapt htseq
module load bowtie2
module load samtools
module load bedtools
module load ucsc_genome

##################  REFERENCE #################
REF=/home/FCAM/aamjad/assembly/Het_chromatin_assembly/dmel_scaffold2_plus0310.fasta

######Index #########
bwa index $REF 

################# REMOVE ADAPTERS AND MAP CUT&Tag DATA WITH BOWTIE2 ##############
for i in Iso1_CIDAM_1
do

### Remove adapter
mkdir fastq_Trim
trim_galore --paired --nextera --length 35 --phred33 --retain_unpaired -o fastq_Trim "$i"_R1.fastq.gz "$i"_R2.fastq.gz

mkdir Mapping
bowtie2 -p 16 -t -x …/index -1 fastq_Trim/"$i"_R1_val_1.fq.gz -2 fastq_Trim/"$i"_R2_val_2.fq.gz | samtools view -bS -@ 16 -f 0x2 | samtools sort -n -m 3G -O bam -o Mapping/"$i"_CENPA-AllReads_Bt2-Dmel_PairsOnly_NameSort.bam
bedtools bamtobed -bedpe -i Mapping/"$i"_CENPA-AllReads_Bt2-Dmel_PairsOnly_NameSort.bam > Mapping/"$i"_CENPA-AllReads_Bt2-Dmel_PairsOnly_NameSort.bed
done

samtools view -@ 16 -h -L ~/assembly/new_annotation/repeat_anno_dmel_new.bed -o Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.bam \
       Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort.bam 
awk '$1==$4 && $6-$2 < 1000 {print $0}' Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.pe.bed > \
       Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.clean.pe.bed
cut -f 1,2,6 Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.clean.pe.bed | sort -k1,1 -k2,2n -k3,3n > \
       Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.frag.pe.bed
bedtools genomecov -bg -i Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.frag.pe.bed -g contigs.size > \
       Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.frag.pe.bedgraph
export LC_COLLATE=C
       sort -k1,1 -k2,2n -k3,3n Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.frag.pe.bedgraph > \
               Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.frag.pe.sort.bedgraph

bedGraphToBigWig Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.frag.pe.sort.bedgraph contigs.size \
       Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.frag.pe.bigwig


#########################################################
#################### Peak Calling #######################
#########################################################
mkdir macs
macs2 callpeak -t Iso1_CIDAM-CENPA-AllReads_bt2-Dmel_PairsOnly_NameSort_Dmel.bam -c CnT_IgG-bt2_PairsOnly_NameSort.bam -f BAMPE -g dm -q 0.01 -B --callsummits

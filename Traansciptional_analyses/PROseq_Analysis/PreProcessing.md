## The following steps were applied for pre-processing and mapping:

1) Run fastqc on raw fastqs (view html file)
2) Trim raw files using CutAdapt (quality cut-off of 20 and length cut-off of 20nt)
3) Run fastqc on trimmed fastqs (view html file) 
4) Get the reverse complement of all clipped reads (ONLY DO THIS FOR PROSEQ) using fastx toolkit
5) Map trimmed and reverse complemented reads to human genome to remove spike-ins using Bowtie2; export as sorted bam using samtools
6) Pull-out only those reads that do NOT map to human genome (aka. Drosophila only reads) using samtools
7) Convert Drosophila only bam into fastq file format using samtools 
8) Map Drosophila-only fastq to Drosophila genome (heterochromatin assembly) using Bowtie2 (with multi-mappers); export as sorted bam using samtools
9) Convert sorted bam to bed file format using bedtools 

## Run pre-processing and mapping code:
Define variables:

    INpath=/path/to/input/directory/
    INfileprefix=SAMPLE1
    DmelanogasterIndex=/labs/Mellone/Asna/Assembly/bowtie2_Index/dmel_scaffold2_plus0310.fasta
    HumangenomeIndex=/isg/shared/databases/alignerIndex/animal/hg38_ucsc/hg38_bowtie2/Human_genome
    OUTpath=/path/to/output/directory/

Run code:

    fastqc ${INpath}${INfileprefix}.fastq
    cutadapt -j 24 -q 20 -m 20 -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC -o ${INpath}${INfileprefix}_R1_trimmed.fastq -p $ {INpath}$   {INfileprefix}_R2_trimmed.fastq ${INpath}${INfileprefix}_R1.fastq ${INpath}${INfileprefix}_R2.fastq
    fastqc ${INpath}${INfileprefix}_trimmed.fastq
    fastx_reverse_complement -Q33 -v -i ${INpath}${INfileprefix}_R1_trimmed.fastq -o ${INpath}${INfileprefix}_R1_trimmed_RC.fastq
    fastx_reverse_complement -Q33 -v -i ${INpath}${INfileprefix}_R2_trimmed.fastq -o ${INpath}${INfileprefix}_R2_trimmed_RC.fastq
    bowtie2 -p 24 -t --very-sensitive -x ${HumangenomeIndex} -1 ${INpath}${INfileprefix}_R1_trimmed_RC.fastq -2 ${INpath}${INfileprefix}_R2_trimmed_RC.fastq | samtools view -@ 24 -bS | samtools sort -O bam -o ${INpath}${INfileprefix}_trimmed_RC_BT2-vs-Hg_sort.bam
    samtools view -@ 20 -b -f 4 -o ${INpath}${INfileprefix}_trimmed_RC_sort_unmapped2HG.bam ${INpath}${INfileprefix}_trimmed_RC_BT2-vs-Hg_sort.bam                  
    samtools fastq ${INpath}${INfileprefix}_trimmed_RC_sort_unmapped.bam > ${INpath}${INfileprefix}_trimmed_RC_sort_unmapped2Hg.fastq
    bowtie2 -p 24 -t -k 100 -x ${DmelanogasterIndex} -1 ${INpath}${INfileprefix}_R1_trimmed_RC_sort_unmapped2Hg.fastq -2 ${INpath}$    {INfileprefix}_R2_trimmed_RC_sort_unmapped2Hg.fastq | samtools view -@ 24 -bS | samtools sort -O bam -o ${OUTpath}${INfileprefix}_trimmed_RC_sort_unmapped2HG_mapped2DM_sort.bam
    bedtools bamtobed -i ${OUTpath}${INfileprefix}_trimmed_RC_sort_unmapped2HG_mapped2DM_sort.bam > ${OUTpath}${INfileprefix}_trimmed_RC_sort_unmapped2HG_mapped2DM_sort.bed


module load bedtools

###### MERYL FOR KMER GENERATION ##########

Local install of Merylv1.3: https://github.com/marbl/meryl
wget https://github.com/marbl/meryl/releases/download/v1.3/meryl-1.3.Linux-amd64.tar.xz
tar -xJf meryl-1.3.Linux-amd64.tar.xz
    
##Generation of unique assembly 51mers:
/home/FCAM/aamjad/meryl-1.3/bin/meryl count k=51 dmel_scaffold2_plus0310.fasta output dmel_het_assembly_51mer_ALLv1.3.meryl 
/home/FCAM/aamjad/meryl-1.3/bin/meryl equal-to 1 dmel_het_assembly_51mer_ALLv1.3.meryl output dmel_het_assembly_51mer_SINGLEv1.3.meryl

##Create unique_51-markers.bed: Compare kmers in input sequences against kmers in input meryl databases and output as a bed. 
/home/FCAM/aamjad/meryl-1.3/bin/meryl-lookup -bed -sequence dmel_scaffold2_plus0310.fasta -mers dmel_het_assembly_51mer_SINGLEv1.3.meryl -output dmel_het_assembly_51mer_SINGLEv1.3.meryl.bed

##Merge overlapping unique 51mers to reduce size of file:
bedtools merge -i dmel_het_assembly_51mer_SINGLEv1.3.meryl.bed > dmel_het_assembly_51mer_SINGLEv1.3.meryl_merge.bed

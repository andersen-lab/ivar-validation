#!/bin/bash

file=$1
prefix=$2

q=15
min_len=30

samtools view -F 4 -h -bS $file "Puerto:360-819" > $prefix.tmp.1.bam
samtools view -F 4 -h -bS $file "Puerto:6674-7096" > $prefix.tmp.2.bam
samtools merge -f $prefix.select_region.bam $prefix.tmp.1.bam $prefix.tmp.2.bam
rm $prefix.tmp.1.bam $prefix.tmp.2.bam

# Fix mate information
samtools sort -n -o $prefix.select_region.tmp.bam $prefix.select_region.bam
samtools fixmate $prefix.select_region.tmp.bam $prefix.select_region.tmp2.bam

# Get only reads mapped in pairs
samtools view -f 1 -h -bS $prefix.select_region.tmp2.bam | samtools sort -n -o $prefix.select_region.name.bam

rm $prefix.select_region.tmp.bam $prefix.select_region.tmp2.bam

java -jar ~/Documents/picard.jar SamToFastq INPUT=$prefix.select_region.name.bam FASTQ=$prefix.select_region.R1.fastq SECOND_END_FASTQ=$prefix.select_region.R2.fastq UNPAIRED_FASTQ=$prefix.select_region.unpaired.fastq

bwa mem ~/Documents/code/ivar/data/db/ZIKV_PRV.fasta $prefix.select_region.R1.fastq $prefix.select_region.R2.fastq | samtools view -F 4 -bS | samtools sort -o $prefix.select_region.realigned.bam

samtools index $prefix.select_region.realigned.bam

# Trimming
ivar trim -b data/test_primer.bed -i $prefix.select_region.realigned.bam -p $prefix.select_region.ivar.trimmed -m $min_len -q $q -s 4 && samtools sort -o $prefix.select_region.ivar.trimmed.sorted.bam $prefix.select_region.ivar.trimmed.bam

cutadapt -e 0.1 -g ^AAGAAAGATCTGGCTGCCATGCT -g ^AGCATGGCAGCCAGATCTTTCTT -g ^ACTAAGGTTGGTCCAAACGCTG -g ^TGATTCCAACCAGGTTTGCGAC -g ^CGTCTTGATGAGGAACAAGGGC -g ^GCCCTTGTTCCTCATCAAGACG -g ^TTCACCAGTGACGTACAACCTG -g ^AAGTGGTCACTGCATGTTGGAC -G ^AAGAAAGATCTGGCTGCCATGCT -G ^AGCATGGCAGCCAGATCTTTCTT -G ^ACTAAGGTTGGTCCAAACGCTG -G ^TGATTCCAACCAGGTTTGCGAC -G ^CGTCTTGATGAGGAACAAGGGC -G ^GCCCTTGTTCCTCATCAAGACG -G ^TTCACCAGTGACGTACAACCTG -G ^AAGTGGTCACTGCATGTTGGAC -m $min_len -q $q,$q -o $prefix.select_region.trimmed.cutadapt.R1.fastq -p $prefix.select_region.trimmed.cutadapt.R2.fastq $prefix.select_region.R1.fastq $prefix.select_region.R2.fastq

# Trim cutadapt
java -jar ~/Documents/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 $prefix.select_region.R1.fastq $prefix.select_region.R2.fastq $prefix.select_region.trimmed.trimmomatic.R1.fastq $prefix.select_region.trimmed.trimmomatic.R1.unpaired.fastq $prefix.select_region.trimmed.trimmomatic.R2.fastq $prefix.select_region.trimmed.trimmomatic.R2.unpaired.fastq  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:$q MINLEN:$min_len HEADCROP:22

# -g ^AAGAAAGATCTGGCTGCCATGCT -g ^AGCATGGCAGCCAGATCTTTCTT -g ^ACTAAGGTTGGTCCAAACGCTG -g ^TGATTCCAACCAGGTTTGCGAC -g ^CGTCTTGATGAGGAACAAGGGC -g ^GCCCTTGTTCCTCATCAAGACG -g ^TTCACCAGTGACGTACAACCTG -g ^AAGTGGTCACTGCATGTTGGAC

# Align cutadapt trimmed
bwa mem ~/Documents/code/ivar/data/db/ZIKV_PRV.fasta $prefix.select_region.trimmed.cutadapt.R1.fastq $prefix.select_region.trimmed.cutadapt.R2.fastq | samtools view -bS -F 4 | samtools sort -o $prefix.select_region.trimmed.cutadapt.sorted.bam

# Align trimmomatic trimmed
bwa mem ~/Documents/code/ivar/data/db/ZIKV_PRV.fasta $prefix.select_region.trimmed.trimmomatic.R1.fastq $prefix.select_region.trimmed.trimmomatic.R2.fastq | samtools view -bS -F 4 | samtools sort -o $prefix.select_region.trimmed.trimmomatic.sorted.bam

samtools depth -d 200000 $prefix.select_region.realigned.bam > $prefix.select_region.depth
samtools depth -d 200000 $prefix.select_region.trimmed.cutadapt.sorted.bam > $prefix.select_region.cutadapt.depth
samtools depth -d 200000 $prefix.select_region.ivar.trimmed.sorted.bam > $prefix.select_region.ivar.depth
samtools depth -d 200000 $prefix.select_region.trimmed.trimmomatic.sorted.bam > $prefix.select_region.trimmomatic.depth

# ./trim_bam.sh data/ZI-merge-26_a.sorted.bam data/ZI-merge-26_a.sorted
# ./trim_bam.sh data/ZI-merge-27_a.sorted.bam data/ZI-merge-27_a.sorted

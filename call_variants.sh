#!/bin/bash

# BC01
time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC01.trimmed.sorted.bam | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2snp --min-coverage 0 --min-reads2 0 --min-avg-qual 0 --min-var-freq 0 > data/variants/BC01.trimmed.sorted.varscan.tsv

time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC01.trimmed.sorted.bam | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2indel --min-coverage 0 --min-reads2 0 --min-avg-qual 0 --min-var-freq 0 > data/variants/BC01.trimmed.sorted.varscan.indel.tsv

time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC01.trimmed.sorted.bam | ivar variants -p data/variants/BC01.trimmed.sorted.ivar -q 0 -t 0

time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC01.trimmed.sorted.bam | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2snp --min-coverage 0 --min-reads2 0 --min-avg-qual 20 min-var-freq 0.03 > data/variants/BC01.trimmed.sorted.varscan.2.tsv

time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC01.trimmed.sorted.bam | ivar variants -p data/variants/BC01.trimmed.sorted.ivar.2.tsv -q 20 -t 0.03


# BC02
time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC02.trimmed.sorted.bam | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2snp --min-coverage 0 --min-reads2 0 --min-avg-qual 0 --min-var-freq 0 > data/variants/BC02.trimmed.sorted.varscan.tsv

time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC02.trimmed.sorted.bam | ivar variants -p data/variants/BC02.trimmed.sorted.ivar -q 0 -t 0

time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC02.trimmed.sorted.bam | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2snp --min-coverage 0 --min-reads2 0 --min-avg-qual 20 min-var-freq 0.03 > data/variants/BC02.trimmed.sorted.varscan.2.tsv
time samtools mpileup --reference data/ref/ZIKV_REF.fasta data/BC02.trimmed.sorted.bam | ivar variants -p data/variants/BC02.trimmed.sorted.ivar.2.tsv -q 20 -t 0.03

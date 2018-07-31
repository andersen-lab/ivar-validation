#!/bin/bash

cd data
wgsim -S 112358 -N 10000 -R 0.2 -r 0.02 -e 0.1 -1 150 -2 150 ref/ZIKV_REF.fasta simulate1.R1.fastq simulate1.R2.fastq
bwa mem ref/ZIKV_REF simulate1.R1.fastq simulate1.R2.fastq | samtools view -bS -F 4 | samtools sort -o simulate1.align.sorted.bam
samtools index simulate1.align.sorted.bam
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate1.align.sorted.bam | ivar consensus -p simulate1.consensus -q 0 -t 0
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate1.align.sorted.bam | ivar consensus -p simulate1.25.consensus -q 0 -t 0.25
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate1.align.sorted.bam | ivar consensus -p simulate1.50.consensus -q 0 -t 0.5
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate1.align.sorted.bam | ivar consensus -p simulate1.90.consensus -q 0 -t 0.90

wgsim -S 112358 -N 10000 -R 0.2 -r 0.01 -e 0.1 -1 150 -2 150 ref/ZIKV_REF.fasta simulate2.R1.fastq simulate2.R2.fastq
bwa mem ref/ZIKV_REF simulate2.R1.fastq simulate2.R2.fastq | samtools view -bS -F 4 | samtools sort -o simulate2.align.sorted.bam
samtools index simulate2.align.sorted.bam
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate2.align.sorted.bam | ivar consensus -p simulate2.consensus -q 0 -t 0
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate2.align.sorted.bam | ivar consensus -p simulate2.25.consensus -q 0 -t 0.25
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate2.align.sorted.bam | ivar consensus -p simulate2.50.consensus -q 0 -t 0.50
samtools mpileup -d 1000 -A -Q 0 -F 0 simulate2.align.sorted.bam | ivar consensus -p simulate2.90.consensus -q 0 -t 0.90

samtools mpileup -d 300000 -A -Q 0 -F 0 ZI-merge-27_a.sorted.bam | ivar consensus -p ZI-merge-27a.consensus -q 0 -t 0
samtools mpileup -d 300000 -A -Q 0 -F 0 ZI-merge-27_a.sorted.bam | ivar consensus -p ZI-merge-27a.25.consensus -q 0 -t 0.25
samtools mpileup -d 300000 -A -Q 0 -F 0 ZI-merge-27_a.sorted.bam | ivar consensus -p ZI-merge-27a.50.consensus -q 0 -t 0.50
samtools mpileup -d 300000 -A -Q 0 -F 0 ZI-merge-27_a.sorted.bam | ivar consensus -p ZI-merge-27a.90.consensus -q 0 -t 0.90

samtools mpileup -d 600000 -A -Q 0 -F 0 ZI-merge-26_a.sorted.bam | ivar consensus -p ZI-merge-26a.consensus -q 0 -t 0
samtools mpileup -d 600000 -A -Q 0 -F 0 ZI-merge-26_a.sorted.bam | ivar consensus -p ZI-merge-26a.25.consensus -q 0 -t 0.25
samtools mpileup -d 600000 -A -Q 0 -F 0 ZI-merge-26_a.sorted.bam | ivar consensus -p ZI-merge-26a.50.consensus -q 0 -t 0.50
samtools mpileup -d 600000 -A -Q 0 -F 0 ZI-merge-26_a.sorted.bam | ivar consensus -p ZI-merge-26a.90.consensus -q 0 -t 0.90

# High Quality for primer trimming
wgsim -S 112358 -N 1000 -1 150 -2 150 ref/ZIKV_REF.fasta simulate3.R1.fastq simulate3.R2.fastq
bwa mem ref/ZIKV_REF simulate3.R1.fastq simulate3.R2.fastq | samtools view -bS -F 4 | samtools sort -o simulate3.align.sorted.bam
samtools index data/simulate3.align.sorted.bam

wgsim -S 112358 -N 1000 -1 150 -2 150 ref/ZIKV_REF.fasta simulate4.R1.fastq simulate4.R2.fastq
bwa mem ref/ZIKV_REF simulate4.R1.fastq simulate4.R2.fastq | samtools view -bS -F 4 | samtools sort -o simulate4.align.sorted.bam

samtools index data/simulate4.align.sorted.bam


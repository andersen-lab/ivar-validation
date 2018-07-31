#!/bin/bash

file=$1				# BAM File
ref=$2				# Reference Fasta
prefix=$3				# Prefix of output

echo "File: " $file
echo "Ref: "$ref
echo "Prefix: "$prefix

# -B is important
samtools mpileup --reference $ref -A -d 600000 -F 0 -B -Q 0 $file > $prefix.pileup

# Qual = 0. Freq = 0.

q=0
f=0

echo "Qual: " $q
echo "Freq: " $f

# VarScan SNP
cat $prefix.pileup | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2snp --min-coverage 0 --min-reads2 0 --min-avg-qual $q --min-var-freq $f > $prefix.varscan.tsv

# VarScan INDEL
cat $prefix.pileup | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2indel --min-coverage 0 --min-reads2 0 --min-avg-qual $q --min-var-freq $f > $prefix.varscan.indel.tsv

# iVar SNP + INDEL
cat $prefix.pileup | ivar variants -p $prefix.ivar -q $q -t $f


# Qual = 20. FREQ 0.03

q=20
f=0.03

echo "Qual: " $q
echo "Freq: " $f

# VarScan SNP
cat $prefix.pileup | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2snp --min-coverage 0 --min-reads2 0 --min-avg-qual $q --min-var-freq $f > $prefix.2.varscan.tsv

# VarScan INDEL
cat $prefix.pileup | java -jar ~/Documents/VarScan.v2.3.9.jar mpileup2indel --min-coverage 0 --min-reads2 0 --min-avg-qual $q --min-var-freq $f > $prefix.2.varscan.indel.tsv

# iVar SNP + INDEL
cat $prefix.pileup | ivar variants -p $prefix.2.ivar -q $q -t $f

rm $prefix.pileup

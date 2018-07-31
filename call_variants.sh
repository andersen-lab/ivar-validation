#!/bin/bash

# 50a
./call_variants_pileup.sh data/ZI-merge-50_a.sorted.bam data/ref/ZIKV_PRV.fasta data/variants/ZI-merge-50_a.sorted

# 49a
./call_variants_pileup.sh data/ZI-merge-49_a.sorted.bam data/ref/ZIKV_PRV.fasta data/variants/ZI-merge-49_a.sorted

# Simuate1
./call_variants_pileup.sh data/simulate1.align.sorted.bam data/ref/ZIKV_REF.fasta data/variants/simulate1.align.sorted

# Simuate2
./call_variants_pileup.sh data/simulate2.align.sorted.bam data/ref/ZIKV_REF.fasta data/variants/simulate2.align.sorted

#!/bin/bash

# BC01
./call_variants_pileup.sh data/BC01.trimmed.sorted.bam data/ref/ZIKV_REF.fasta data/variants/BC01.trimmed.sorted

# BC02
./call_variants_pileup.sh data/BC02.trimmed.sorted.bam data/ref/ZIKV_REF.fasta data/variants/BC02.trimmed.sorted

# 26a
./call_variants_pileup.sh data/ZI-merge-26_a.sorted.bam data/ref/ZIKV_PRV.fasta data/variants/ZI-merge-26_a.sorted

# 27a
./call_variants_pileup.sh data/ZI-merge-27_a.sorted.bam data/ref/ZIKV_PRV.fasta data/variants/ZI-merge-27_a.sorted

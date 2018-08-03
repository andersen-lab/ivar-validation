# Validation of iVar

Repository contains scripts used to validation of trimming, consensus calling and variant calling using [iVar](https://github.com/andersen-lab/ivar/). 

```
.
├── data
│   ├── alignments - Contains alignments of consensus sequences called using iVar, Geneious
│   ├── ref - Reference sequences used
│   ├── trimmomatic - Trimmed fastq files using Trimmomatic
│   └── variants - Variant csv and tsv files called using VarScan and iVar.
├── plots - Plots generated by R scripts and used to create validation figures.
└── scripts - Contains R scripts used to create validation figures and validate results. 
```

### Details:

| Command | Validated Against |
|:--------|:------------------|
|  trim   | [cutadapt](https://github.com/marcelm/cutadapt/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)|
| consensus | [Geneious](https://www.geneious.com/) |
| variants | [VarScan2](https://dkoboldt.github.io/varscan/) |
# SeqAna_Pipe

A Nextflow-based RNA-seq analysis pipeline for processing raw FASTQ files and generating gene count matrices, with downstream differential expression analysis for centriolar mRNA localization studies.

## Overview

This pipeline repo includes RNA-seq data processing from raw reads to gene counts, with specialized analysis workflows for identifying centriole-enriched transcripts. It includes quality control, alignment, quantification, and downstream analysis notebooks.

## Pipeline details

- **Quality Control**: FastQC analysis of raw reads
- **Alignment**: STAR aligner with automatic index generation
- **Quantification**: featureCounts for gene-level counting
- **Reporting**: MultiQC aggregated quality reports
- **Downstream Analysis**: Jupyter notebooks for differential expression and visualization

## Dependencies
- Nextflow ≥ 23.04.0
- FastQC 0.12.1
- STAR 2.7.10b
- Subread 2.0.6 (featureCounts)
- SAMtools 1.17
- MultiQC 1.15
- Python 3.9
## Packgaes for downstream analysis
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn
- mygene (annotate genes)

## Run the pipeline
Install via conda/mamba using the provided environment file:

### prepare the sample sheet
```bash
sample_id,read1,read2,cell_line,replicate
rep1_0,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,HeLa,1
rep1_1,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,HeLa,2
```

### Run
```bash
conda env create -f environment.yml
conda activate rnaseq-pipeline

nextflow run main.nf \
  --samplesheet samples.csv \
  --genome reference.fa \
  --gtf annotations.gtf \
  --outdir results
```


## Downstream Analysis

### Analysis Scripts

1. RNAseq_Analysis.ipynb: Main differential expression analysis

    - Identifies centriole-enriched mRNAs
    - Removes cell line-specific genes
    - PCA visualization
    - MA plots for differential expression

2. variation_between_replicates.ipynb: Quality assessment
    - Replicate correlation analysis
    - Expression distribution comparisons
    - Statistical validation

3. analysis.py: Automated preprocessing script

4. mapGeneID.ipynb: map Gene symbol annotation using MyGene


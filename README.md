# ChIPseq_humanrDNA

This repository contains ChIP-seq data analysis pipeline for ribosomal DNA (rDNA) repeats and genome-wide (For paired-end short read - Illumina data). 

A custom manually masked reference genome was utilised in this analysis to account for the highly repetitive nature of rDNA. 
The workflow followed in creating a custom reference is outlined in https://github.com/tudumanne/custom-reference-human.

### Pipeline overview

The pipeline was run on an HPC (high-performance computing) system based on CentOS (Linux). The 'scripts' folder contains template bash scripts.

### Software installation 

The required command-line tools were installed via conda on Linux. 

Miniconda documentation https://docs.conda.io/en/latest/miniconda.html

- Create an environment named chip-seq using a .yaml file, which defines the tools that need to be installed. 

```console
conda env create -n chip-seq -f environment.yaml
```
  
1. Quality check of raw fastq files - FastQC/MultiQC

```console
fastqc -o fastqc --extract --dir fastqc_chip --format fastq h3k4me3_*.fastq.gz
fastqc -o fastqc --extract --dir fastqc_input --format fastq input_*.fastq.gz

multiqc fastqc_chip/
```

2. Read alignment, processing and post-alignment quality check - Bowtie2 and Samtools 

Read alignment using Bowtie2, SAM to BAM conversion, sorting and indexing the BAM files using Samtools

```console
bowtie2 -x index -1 {sample}_R1.fastq.gz -2 {sample}_R2.fastq.gz | samtools view -bS > {sample}.bam
samtools sort {sample}.bam > {sample}_sorted.bam
samtools index {sample}_sorted.bam
```

Quality check of BAM files

```console
samtools stats
```

3. Coverage track generation and visualisation - deepTools and IGV

Generate a coverage track 

```console
bamCoverage -b {sample}_sorted.bam -o {sample}_coverage.bw --outFileFormat bigwig --binSize 10 --smoothLength 25 --extendReads --ignoreForNormalization chrX chrY chrMT 
```

Generate a coverage track normalised to input

```console
bamCompare -b1 {sample}_sorted.bam -b2 input_merged.bam -o {sample}.bw --scaleFactorsMethod None --operation ratio --binSize 10 --normalizeUsing RPKM --smoothLength 25 --extendReads --ignoreForNormalization chrX chrY chrMT
```

4. Peak calling - MACS2

```console
macs2 {sample}_sorted.bam -c {input_merged}.bam -g hs -n sample -q 0.05 --cutoff -f BAMPE --keep-dup all
```

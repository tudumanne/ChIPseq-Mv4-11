# ChIPseq_humanrDNA

This repository contains ChIP-seq data analysis pipeline for ribosomal DNA (rDNA) repeats and genome-wide (For paired-end short read - Illumina data). 

A custom manually masked reference genome was utilised in this analysis to account for the highly repetitive nature of rDNA. 

# custom-reference-human
Custom human reference genome for ribosomal DNA (rDNA) analysis

### Workflow - creating a custom manually masked reference 
<img width="500" alt="image" src="https://github.com/user-attachments/assets/1e89e246-6e82-42af-9840-0200d4cfa6ee" />

#### Identify and mask regions across the human genome (GRCh38) that can interfere with read mapping across rDNA during ChIP-seq data analysis

1. First the rDNA canonical copy (KY962518.1) was fragmented in-silico to generate 150bp long overlapping fragments with step size of 1bp.

  - 150bp fragment length was selected considering the average fragment size of input samples which is between 200-300bp. 
  - Example tool that can be used for in-silico fragment generation - https://www.genecorner.ugent.be/split_fasta.html  	

2. Generated fragments were aligned to the standard unmasked reference (without incorporating rDNA canonical sequence) using Bowtie2, allowing up to 1000 possible alignments (k=1000)

- Considering the possible sequence variation in experimental samples stringent settings (alignment without rDNA copy and allowing multiple alignments) were used to identify ‘rDNA-like’ regions.

```console
#build index
bowtie2-build GRCh38_unmasked_chr_rdnamasked_rdna_edited.fa GRCh38_unmasked_chr_rdnamasked_rdna_edited_index -p 24

#align the in-silico fragmented reads to the reference and convert to BAM format
bowtie2 -x GRCh38_unmasked_chr_rdnamasked_rdna_edited_index -U rDNA_150bp.fa -f -p 24 -k 1000 | samtools view -bS > unmasked_150bp.bam

#sort the BAM file and index
samtools sort unmasked_150bp.bam > unmasked_150bp_sorted.bam
samtools index unmasked_150bp_sorted.bam
```

3. ‘rDNA-like’ sequences across the genome were extracted as a .bed file (scripts folder - circos_plot.R) and the identified regions were masked in the GRCh38 reference (without scaffolds) using 'maskfasta' - bedtools.

```console
bedtools maskfasta -fi reference.fa -bed regions.bed -fo reference_rdna_masked.fa
```

Mappability percentage and mapping quality (MAPQ) across rDNA for different reference genomes were determined using Samtools 'mpileup' command.

```console
samtools mpileup -f reference.fa -B -r rDNA_repeat -o input_mpileup -O -s -a -R input_seq.bam
```

- Canonical rDNA reference was incorporated to the resulting fasta file which was then used as the custom reference for ChIP-seq data analysis.

### ChIP-seq data analysis pipeline overview

The pipeline was run on an HPC (high-performance computing) system based on CentOS (Linux).

### Software installation 

The required command-line tools were installed via conda on Linux. 

Miniconda documentation https://docs.conda.io/en/latest/miniconda.html

- Create an environment named chip-seq using a .yaml file, which defines the tools that need to be installed. 

```console
conda env create -n chip-seq -f environment.yaml
```
  
1. Quality check of raw fastq files - FastQC/MultiQC

```console
#run fastqc 
for i in *_001.fastq.gz;
  do name=$(basename ${i} _001.fastq.gz);
  fastqc -o fastqc --extract --dir fastqc_reports --format fastq ${name}_001.fastq.gz;
done

#compile a report using multiqc
multiqc fastqc_reports/
```

2. Read alignment, processing and post-alignment quality check - Bowtie2 and Samtools 

Read alignment using Bowtie2, SAM to BAM conversion, sorting and indexing the BAM files using Samtools

```console
#build index
bowtie2-build GRCh38_unmasked_chr_rdnamasked_rdna_edited.fa GRCh38_unmasked_chr_rdnamasked_rdna_edited_index -p 24

#perform alignment
bowtie2 -x GRCh38_unmasked_chr_rdnamasked_rdna_edited_index -1 {sample}_R1_001.fastq.gz -2 {sample}_R2_001.fastq.gz -p 24 -S {sample}.sam

#convert to BAM
samtools view -S -b {sample}.sam > {sample}.bam

#sort and index
samtools sort -T {sample}.bam > {sample}_sorted.bam
samtools index {sample}_sorted.bam

```

Quality check of BAM files

```console
samtools stats {sample}_sorted.bam
```

Merge input BAM files and index

```console
samtools merge input1_merged.bam input1_n1_sorted.bam input1_n2_sorted.bam input1_n3_sorted.bam
samtools index /g/data/ff84/NH_ChIPseq_Feb2022/fastq_files/input3_n3_sorted.bam
```

3. Coverage track generation and visualisation - deepTools and IGV

Generate a coverage track normalised to input

```console
bamCompare -b1 {sample}_sorted.bam -b2 input_merged.bam -o {sample}.bw --scaleFactorsMethod None --outFileFormat bigwig --binSize 10 --smoothLength 25 --extendReads --normalizeUsing RPKM --ignoreForNormalization chrX chrY KY962518_edited chrMT -p 24 --operation ratio
```

4. Peak calling - MACS2

```console
macs2 {sample}_sorted.bam -c {input_merged}.bam -g hs -n sample -q 0.05 --cutoff -f BAMPE --keep-dup all
```

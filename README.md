# ChIPseq-Mv4-11 code repository

This repository contains ChIP-seq data analysis pipeline for ribosomal DNA (rDNA) repeats and genome-wide for paired-end short read (Illumina) data. 

Note: A custom manually masked reference genome was utilised in this analysis to account for the highly repetitive nature of rDNA. 

## Creation of a custom human reference genome for ribosomal DNA (rDNA) analysis

Due to its highly repetitive nature rDNA loci are not included in the human genome assemblies to date and two canonical reference sequences are available through
GenBank accession no.  KY962518 and U13369. There are several rDNA pseudogenes or ‘rDNA-like’ fragments present across the in-silico human reference genome which interfere with short read alignment. Therefore, rDNA chromatin structure characterisation using ChIP-seq requires the development and testing of novel bioinformatics analysis methods to ensure robust quantification of specific enrichment at rDNA loci.

### Workflow - creating a custom manually masked reference 
<img width="500" alt="image" src="https://github.com/user-attachments/assets/1e89e246-6e82-42af-9840-0200d4cfa6ee" />

#### Identify and mask regions across the human genome (GRCh38) that can interfere with read mapping across rDNA during ChIP-seq data analysis

1. First the rDNA canonical copy (KY962518.1) was fragmented in-silico to generate 150bp long overlapping fragments with step size of 1bp.

  - 150bp fragment length was selected considering the average fragment size of input samples, which is between 200-300bp. 
  - Example tool that can be used for in-silico fragment generation - https://www.genecorner.ugent.be/split_fasta.html  	

2. Generated fragments were aligned to the standard unmasked reference (without incorporating rDNA canonical sequence) using Bowtie2 v2.3.5.1, allowing up to 1000 possible alignments (k=1000)

- Considering the possible sequence variation in experimental samples, stringent settings (alignment without rDNA copy and allowing multiple alignments) were used to identify ‘rDNA-like’ regions.

```console
#build index
bowtie2-build unmasked_grch38.fa unmasked_grch38_index

#align the in-silico fragmented reads to the reference and convert to BAM format using Samtools v2.30.0
bowtie2 -x unmasked_grch38_index -U rDNA_150bp.fa -f -k 1000 | samtools view -bS > unmasked_150bp.bam

#sort the BAM file and index
samtools sort unmasked_150bp.bam > unmasked_150bp_sorted.bam
samtools index unmasked_150bp_sorted.bam
```

3. ‘rDNA-like’ sequences mapped across the genome were extracted as a .bed file using R 4.1 (scripts folder - extract_genomic_regions.R).
   
5.  The identified regions were masked in the GRCh38 reference (without scaffolds) using 'maskfasta' tool in BEDTools suite version 2.30.0.

```console
bedtools maskfasta -fi unmasked_grch38.fa -bed regions.bed -fo reference_rdna_masked.fa
```

Mappability percentage and mapping quality (MAPQ) across rDNA for different reference genomes (unmasked vs. hard-masked vs. manually masked reference created) were determined using Samtools (v2.30.0) 'mpileup' command.

```console
samtools mpileup -f reference_rdna_masked.fa -B -r rDNA_repeat -o input_mpileup -O -s -a -R input_seq.bam
```

- Canonical rDNA reference was incorporated to the resulting fasta file which was then used as the custom reference for ChIP-seq data analysis.


### ChIP-seq data analysis pipeline overview

The pipeline was run on an HPC (high-performance computing) system based on CentOS (Linux).

1. Quality check of raw fastq files - FastQC/MultiQC

```console
#!/bin/bash

#run fastqc
for i in *_001.fastq.gz; do
  name=$(basename ${i} _001.fastq.gz)
  fastqc -o fastqc --extract --dir fastqc_reports --format fastq ${name}_001.fastq.gz
done

#compile a report using multiqc
multiqc fastqc_reports/
```

2. Read alignment, processing and post-alignment quality check - Bowtie2 v2.3.5.1 and Samtools v2.30.0

Read alignment using Bowtie2, followed by SAM to BAM conversion, sorting and indexing the BAM files using Samtools

```console
#!/bin/bash

#build bowtie2 index

bowtie2-build reference_rdna_masked.fa reference_rdna_masked_index

#perform alignment

INPUT_DIR="/path/to/fastq"
OUTPUT_DIR="/path/to/output"

for R1 in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="$INPUT_DIR/${SAMPLE}_R2_001.fastq.gz"
    bowtie2 -x reference_rdna_masked_index -1 "$R1" -2 "$R2" | samtools view -bS - > "$OUTPUT_DIR/${SAMPLE}.bam"
done

#sort and index

for bam in *.bam; do
    base=$(basename "$bam" .bam)
    sorted="sorted_bam/${base}.sorted.bam"
    samtools sort -T "$bam" > "$sorted"
    samtools index "$sorted"
done
```

Quality check of BAM files

```console
#!/bin/bash

for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    samtools stats "$bam" > "bam_stats/${sample}.stats"
done
```

Merge input BAM files, sort and index

```console
#!/bin/bash

for INPUT in input1 input2 input3; do
    REPLICATES=$(ls ${INPUT}_rep*.bam | tr '\n' ' ')
    OUT="merged_inputs/${INPUT}_merged.bam"
    samtools merge "$OUT" $REPLICATES
    samtools sort -T "$OUT" > "merged_controls/${INPUT}_merged.sorted.bam" 
    samtools index "merged_controls/${INPUT}_merged.sorted.bam"
done
```

3. Coverage track generation and visualisation - deepTools and IGV

Generate a coverage track normalised to merged input (n=3) 

```console
#!/bin/bash

CHIP_DIR="/path/to/inputsamples"
INPUT_DIR="/path/to/chipedsamples"
OUTPUT_DIR="/path/to/normalised_coverage_output"

for INPUT in "$INPUT_DIR"/input*.bam; do
    INPUT_BASENAME=$(basename "$INPUT" .bam)

  # Run bamCompare

    for i in {1..3}; do
        CHIP="${CHIP_DIR}/${INPUT_BASENAME}_chip${i}.bam"
        OUTFILE="${OUTPUT_DIR}/${INPUT_BASENAME}_chip${i}_ratio.bw"

        bamCompare -b1 "$CHIP" -b2 "$INPUT" \
            --operation ratio \
            -o "$OUTFILE" \
            --outFileFormat bigwig --binSize 10 --smoothLength 25 --normalizeUsing RPKM --extendReads \
            --ignoreForNormalization chrX chrY KY962518_edited chrMT \
            --scaleFactorsMethod None
    done
done
```

4. Peak calling - MACS2

```console
#!/bin/bash

CHIP_DIR="/path/to/inputsamples"
INPUT_DIR="/path/to/chipedsamples"
OUTPUT_DIR="/path/to/peak_calling_output"

for INPUT in "$INPUT_DIR"/input*.bam; do
    INPUT_NAME=$(basename "$INPUT" .bam)

  # Run MACS2 peak calling
    for i in {1..3}; do
        CHIP="${CHIP_DIR}/${INPUT_NAME}_chip${i}.bam"
        OUT_PREFIX="${INPUT_NAME}_chip${i}"
        OUT_DIR="${OUTPUT_DIR}/${OUT_PREFIX}"

        macs2 callpeak -t "$CHIP" -c "$INPUT" \
            -f BAMPE -g hs -n "$OUT_PREFIX" \
            --outdir "$OUT_DIR" \
            -q 0.05 --keep-dup all
    done
done
```

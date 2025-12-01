# GATK Variant Calling Pipeline

This repository contains a complete workflow for performing variant calling using the Genome Analysis Toolkit (GATK).
The pipeline covers all major steps from raw FASTQ files to final annotated variants (VCF).

# Overview

This pipeline demonstrates how to perform:

>Quality control

>Read trimming

>Genome alignment

>BAM processing (sorting, marking duplicates, indexing)

>Base Quality Score Recalibration (BQSR)

>Variant calling (VCF / GVCF)

>Joint genotyping (for multi-sample projects)

>Variant filtering (VQSR / hard filtering)

>Variant annotation (Funcotator / VEP / ANNOVAR)

The instructions work for human datasets (hg38) and can be adapted for bacterial/fungal genomes.

# Environment Setup

Create and activate a dedicated conda environment:

conda create -n gatk_env python=3.9

conda activate gatk_env

conda install -c bioconda gatk4

Verify installation:

gatk --version

# Download Reference Genome
wget -P reference/hg38/ http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip reference/hg38/hg38.fa.gz

# Quality Check — FastQC
fastqc sample_1.fastq sample_2.fastq

# Read Trimming — Trimmomatic
java -jar trimmomatic-0.30.jar PE \
sample_1.fastq sample_2.fastq \
sample_1_paired.fq.gz sample_1_unpaired.fq.gz \
sample_2_paired.fq.gz sample_2_unpaired.fq.gz \
ILLUMINACLIP:adapter.fa:2:30:10 LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:15 MINLEN:36

# Alignment — BWA

Install BWA:

sudo apt install bwa

Index reference:

bwa index reference/hg38/reference_genome.fasta

# Align:

bwa mem reference/hg38/reference_genome.fasta \
sample_1.fastq sample_2.fastq > aligned_reads.sam

# Convert and sort:

samtools view -bS aligned_reads.sam > aligned_reads.bam

java -jar picard.jar SortSam -INPUT aligned_reads.bam -OUTPUT sorted_reads.bam -SORT_ORDER coordinate
samtools index sorted_reads.bam

# Mark Duplicates — Picard
java -jar picard.jar MarkDuplicates \
INPUT=with_readgroups.bam OUTPUT=marked_reads.bam METRICS_FILE=dup_metrics.txt


If required:

java -jar picard.jar AddOrReplaceReadGroups \
-I sorted_reads.bam -O with_readgroups.bam \
-RGID 1 -RGLB library1 -RGPL illumina -RGPU unit1 -RGSM sample1

# Base Quality Score Recalibration (BQSR)

Index reference:

samtools faidx reference_genome.fasta

gatk CreateSequenceDictionary -R reference_genome.fasta


# Download known sites:

wget -P reference/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

wget -P reference/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


# Run BQSR:

gatk BaseRecalibrator \
-I marked_reads.bam -R reference_genome.fasta \
-known-sites dbsnp.vcf -O recal_data.table

gatk ApplyBQSR \
-I marked_reads.bam -R reference_genome.fasta \
--bqsr-recal-file recal_data.table -O recalibrated_reads.bam

# Variant Calling — HaplotypeCaller
Standard VCF:
gatk HaplotypeCaller \
-R reference_genome.fasta \
-I recalibrated_reads.bam \
-O raw_variants.vcf

# GVCF mode:
gatk HaplotypeCaller \
-R reference_genome.fasta \
-I recalibrated_reads.bam \
-O sample.g.vcf -ERC GVCF

# Joint Genotyping (Multi-sample)
gatk CombineGVCFs \
-R reference_genome.fasta -V sample1.g.vcf -V sample2.g.vcf \
-O combined.g.vcf

gatk GenotypeGVCFs \
-R reference_genome.fasta -V combined.g.vcf \
-O combined_variants.vcf

# Variant Filtering
Human data — VQSR:
gatk VariantRecalibrator ...
gatk ApplyVQSR ...

# Bacterial/Fungal — Hard Filtering:
gatk VariantFiltration \
-R reference.fasta -V raw_variants.vcf \
--filter-name "LowQual" \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
-O filtered_variants.vcf


# Separate SNP / INDEL:

gatk SelectVariants ...

# Variant Annotation
Funcotator:
gatk Funcotator \
-R reference/hg38/reference_genome.fasta \
-V filtered_variants.vcf \
--output-file annotated_variants.vcf \
--data-sources-path /path/to/funcotator_data_sources

# VEP (Ensembl)

Alternatively, annotation can be performed through Ensembl VEP (web or CLI).

# For your queries and suggestions

Reach out to 'sajeevrajssr@gmail.com'

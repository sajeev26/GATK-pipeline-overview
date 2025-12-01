#!/bin/bash
# ================================================================
#      GATK Variant Calling Pipeline — Whole Genome / Exome
# ================================================================
# Author : Sajeev Raj S  
# Description : Fully automated pipeline for FASTQ → VCF  
# Dependencies : GATK4, Picard, BWA, SAMtools, Trimmomatic, FastQC
# ================================================================

### -----------------------------
### 1. Create Conda Environment
### -----------------------------
echo "Activating GATK environment..."
conda create -n gatk_env python=3.9 -y
source activate gatk_env
conda install -c bioconda gatk4 fastqc samtools bwa trimmomatic -y

### -----------------------------
### 2. Download Reference Genome
### -----------------------------
mkdir -p reference/hg38
wget -P reference/hg38/ http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip reference/hg38/hg38.fa.gz
REFERENCE="reference/hg38/hg38.fa"

### -----------------------------
### 3. FastQC Quality Control
### -----------------------------
echo "Running FastQC..."
mkdir -p fastqc_output
fastqc sample_1.fastq sample_2.fastq -o fastqc_output/

### -----------------------------
### 4. Read Trimming (Trimmomatic)
### -----------------------------
echo "Trimming reads..."
java -jar trimmomatic-0.30.jar PE \
sample_1.fastq sample_2.fastq \
sample_1_paired.fq.gz sample_1_unpaired.fq.gz \
sample_2_paired.fq.gz sample_2_unpaired.fq.gz \
ILLUMINACLIP:adapter.fa:2:30:10 LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:15 MINLEN:36

### -----------------------------
### 5. Alignment using BWA
### -----------------------------
echo "Indexing reference genome..."
bwa index $REFERENCE

echo "Aligning reads..."
bwa mem $REFERENCE \
sample_1_paired.fq.gz sample_2_paired.fq.gz > aligned_reads.sam

### -----------------------------
### 6. SAM → BAM conversion
### -----------------------------
samtools view -bS aligned_reads.sam > aligned_reads.bam

### -----------------------------
### 7. Sort BAM using Picard
### -----------------------------
echo "Sorting BAM..."
java -jar picard.jar SortSam \
INPUT=aligned_reads.bam \
OUTPUT=sorted_reads.bam \
SORT_ORDER=coordinate

samtools index sorted_reads.bam

### -----------------------------
### 8. Add Read Groups (Required for GATK)
### -----------------------------
java -jar picard.jar AddOrReplaceReadGroups \
-I sorted_reads.bam \
-O with_readgroups.bam \
-RGID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM sample1

samtools index with_readgroups.bam

### -----------------------------
### 9. Mark Duplicates
### -----------------------------
echo "Marking duplicates..."
java -jar picard.jar MarkDuplicates \
INPUT=with_readgroups.bam \
OUTPUT=marked_reads.bam \
METRICS_FILE=dup_metrics.txt

samtools index marked_reads.bam

### -----------------------------
### 10. Base Quality Score Recalibration
### -----------------------------
echo "Generating FASTA index and dictionary..."
samtools faidx $REFERENCE
gatk CreateSequenceDictionary -R $REFERENCE

echo "Downloading known sites..."
wget -P reference/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P reference/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

DBSNP="reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"

### Base Recalibration
gatk BaseRecalibrator \
-R $REFERENCE \
-I marked_reads.bam \
-known-sites $DBSNP \
-O recal_data.table

gatk ApplyBQSR \
-R $REFERENCE \
-I marked_reads.bam \
--bqsr-recal-file recal_data.table \
-O recalibrated_reads.bam

### -----------------------------
### 11. Variant Calling
### -----------------------------
echo "Running HaplotypeCaller..."
gatk HaplotypeCaller \
-R $REFERENCE \
-I recalibrated_reads.bam \
-O raw_variants.vcf

### -----------------------------
### 12. GVCF mode (optional)
### -----------------------------
gatk HaplotypeCaller \
-R $REFERENCE \
-I recalibrated_reads.bam \
-O sample.g.vcf \
-ERC GVCF

### -----------------------------
### 13. Variant Filtering (Hard Filter)
### -----------------------------
gatk VariantFiltration \
-R $REFERENCE \
-V raw_variants.vcf \
--filter-name "LowQual" \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
-O filtered_variants.vcf

### -----------------------------
### 14. Extract SNPs and INDELs
### -----------------------------
gatk SelectVariants \
-R $REFERENCE \
-V filtered_variants.vcf \
--select-type SNP \
-O snps.vcf

gatk SelectVariants \
-R $REFERENCE \
-V filtered_variants.vcf \
--select-type INDEL \
-O indels.vcf

### -----------------------------
### 15. Variant Annotation (Funcotator)
### -----------------------------
# NOTE: User must provide funcotator data sources.
echo "Annotation step requires Funcotator data sources."

# Example:
# gatk Funcotator \
# -R $REFERENCE \
# -V filtered_variants.vcf \
# --output-file annotated_variants.vcf \
# --data-sources-path funcotator_datasources \
# --ref-version hg38

echo "Pipeline completed successfully!"

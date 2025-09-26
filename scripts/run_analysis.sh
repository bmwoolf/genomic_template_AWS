#!/bin/bash
# WGS analysis pipeline
# step 1: install tools

set -e

# install conda if not present
if ! command -v conda &> /dev/null; then
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    rm Miniconda3-latest-Linux-x86_64.sh
fi

# configure conda channels
conda config --add channels bioconda
conda config --add channels conda-forge

# install bioinformatics tools
conda install -y htslib bwa samtools bedtools gatk4 manta

# install Python packages
conda install -y pandas numpy matplotlib seaborn jupyter pysam biopython

# verify installations
htsfile --version
bwa version
samtools --version
bedtools --version
gatk --version
python -c "import pandas, numpy, matplotlib, seaborn, pysam, Bio; print('Python packages OK')"

# step 2: align reads
# download reference genome
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipeline.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# index reference genome
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# align FASTQ to BAM
NCPUS=16
ID=HM1A2ABC1:4
PU=HM1A2ABC1:4:MY_SAMPLE_1
LB=20200101
PL=ILLUMINA
SM=MY_SAMPLE_1
read_group_id="@RG\tID:$ID\tPU:$PU\tLB:$LB\tPL:$PL\tSM:$SM"

bwa mem -t $NCPUS -R "$read_group_id" GCA_000001405.15_GRCh38_no_alt_analysis_set.fna my_fastq.R1.fastq my_fastq.R2.fastq | \
    samtools view -1 - -o my_genome.bam

# step 3: improve quality
# mark duplicates
gatk MarkDuplicates \
    --INPUT my_genome.bam \
    --OUTPUT my_genome_markdup.bam \
    --METRICS_FILE markdup_metrics.txt \
    --VALIDATION_STRINGENCY SILENT

# sort by coordinate
gatk SortSam \
    --INPUT my_genome_markdup.bam \
    --OUTPUT my_genome_sorted.bam \
    --SORT_ORDER coordinate

# index BAM
samtools index my_genome_sorted.bam

# download BQSR resources
wget -q https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
wget -q https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
wget -q https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget -q https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
wget -q https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -q https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# build recalibration model
gatk BaseRecalibrator \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -I my_genome_sorted.bam \
    --known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known-sites Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O bqsr_table.table

# apply recalibration
gatk ApplyBQSR \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -I my_genome_sorted.bam \
    --bqsr-recal-file bqsr_table.table \
    -O my_genome_bqsr.bam

# step 4: call variants
# create gVCF
gatk HaplotypeCaller \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -I my_genome_bqsr.bam \
    -D Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -ERC GVCF \
    -O my_genome.g.vcf

# genotype gVCF
gatk GenotypeGVCFs \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -V my_genome.g.vcf \
    -D Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -O my_variants.vcf

# step 5: annotate results
# install VEP
conda install -y ensembl-vep

# run VEP annotation
vep --input_file my_variants.vcf \
    --output_file my_variants_annotated.vcf \
    --format vcf \
    --vcf \
    --symbol \
    --terms SO \
    --hgvs \
    --fasta GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --offline \
    --cache \
    --merged

echo "~analysis complete: my_variants_annotated.vcf"
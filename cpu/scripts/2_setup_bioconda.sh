#!/bin/bash

# =============================================================================
# BIOCONDA ENVIRONMENT SETUP
# sets up conda environment and downloads reference data
# run this once to prepare the environment for genomic analysis
# takes about 1.5 hours for all downloading and indexing
# =============================================================================

set -euo pipefail

# configuration
REFERENCE_GENOME="GRCh38"
WORK_DIR="${WORK_DIR:-$HOME/wgs_data}"
REF_DIR="$WORK_DIR/references"

# colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # no color

# logging function
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# setup conda environment
setup_environment() {
    log "Setting up conda environment..."
    
    # install miniconda if not present
    if ! command -v conda &> /dev/null; then
        log "Installing miniconda..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
        export PATH="$HOME/miniconda3/bin:$PATH"
        echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    fi
    
    # check if genomics environment exists, create if not
    if ! conda env list | grep -q "genomics"; then
        log "Creating genomics conda environment..."
        conda create -n genomics python=3.9 -y
    else
        info "Genomics environment already exists"
    fi
    
    # initialize conda for this shell and activate environment
    eval "$(conda shell.bash hook)"
    conda activate genomics
    
    # configure conda channels for better compatibility
    log "Configuring conda channels..."
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --set channel_priority strict
    
    # install bioinformatics tools (excluding multiqc due to Python 3.9 conflicts)
    log "Installing compiled bioinformatics tools via conda..."
    conda install -y \
        bwa \
        samtools \
        bedtools \
        gatk4 \
        bcftools \
        htslib \
        picard \
        fastqc \
        trimmomatic \
        pigz \
        pbzip2 \
        ncurses
    
    # install multiqc via pip for better Python 3.9+ compatibility
    log "Installing multiqc via pip..."
    pip install multiqc
    
    # fix ncurses library compatibility for older samtools
    log "Fixing ncurses library compatibility..."
    ln -sf $CONDA_PREFIX/lib/libncurses.so.6 $CONDA_PREFIX/lib/libncurses.so.5 2>/dev/null || true
    ln -sf $CONDA_PREFIX/lib/libtinfo.so.6 $CONDA_PREFIX/lib/libtinfo.so.5 2>/dev/null || true
    
    log "Environment setup complete"
}

# download reference genome
download_reference() {
    log "Downloading reference genome..."
    
    mkdir -p "$REF_DIR"
    cd "$REF_DIR"
    
    # download GRCh38 reference genome
    if [ ! -f "GRCh38.fna" ]; then
        log "Downloading GRCh38 reference genome..."
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip hg38.fa.gz
        mv hg38.fa GRCh38.fna
    else
        info "Reference genome already downloaded"
    fi
    
    # index reference genome
    if [ ! -f "GRCh38.fna.bwt" ]; then
        log "Indexing reference genome with BWA..."
        bwa index GRCh38.fna
    else
        info "BWA index already exists"
    fi
    
    # create FASTA index
    if [ ! -f "GRCh38.fna.fai" ]; then
        log "Creating FASTA index..."
        samtools faidx GRCh38.fna
    else
        info "FASTA index already exists"
    fi
    
    # create sequence dictionary
    if [ ! -f "GRCh38.dict" ]; then
        log "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R GRCh38.fna -O GRCh38.dict
    else
        info "Sequence dictionary already exists"
    fi
    
    log "Reference genome setup complete"
}

# download GATK resources
download_gatk_resources() {
    log "Downloading GATK resources..."
    
    mkdir -p "$REF_DIR/gatk_resources"
    cd "$REF_DIR/gatk_resources"
    
    # download BQSR resources
    if [ ! -f "Homo_sapiens_assembly38.known_indels.vcf.gz" ]; then
        log "Downloading known indels..."
        wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
        wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
    else
        info "Known indels already downloaded"
    fi
    
    if [ ! -f "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ]; then
        log "Downloading Mills indels..."
        wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
    else
        info "Mills indels already downloaded"
    fi
    
    if [ ! -f "Homo_sapiens_assembly38.dbsnp138.vcf.gz" ]; then
        log "Downloading dbSNP..."
        wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
        wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
    else
        info "dbSNP already downloaded"
    fi
    
    log "GATK resources download complete"
}

# main execution
main() {
    log "==================================================================="
    log "Starting bioconda environment setup..."
    log "==================================================================="
    log "Reference: $REFERENCE_GENOME"
    log "Work directory: $WORK_DIR"
    
    # create directory structure
    mkdir -p "$WORK_DIR" "$REF_DIR"
    
    # run setup steps
    setup_environment
    download_reference
    download_gatk_resources
    
    log "==================================================================="
    log "Setup complete! Environment ready for analysis."
    log "==================================================================="
    log "Next step: Run './scripts/workstation/analyze_genome.sh' to start analysis."
}

# run main function
main "$@"
    
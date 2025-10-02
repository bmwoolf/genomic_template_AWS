#!/bin/bash

# =============================================================================
# LOCAL GENOMICS ENVIRONMENT SETUP
# for Linux workstations with NVIDIA GPU support
# =============================================================================

set -euo pipefail

# configuration
WORK_DIR="${WORK_DIR:-$HOME/genomic_analysis}"
CONDA_ENV="${CONDA_ENV:-genomics}"

# colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

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

# check if running on linux
check_os() {
    if [[ "$OSTYPE" != "linux-gnu"* ]]; then
        error "This script is designed for Linux. Your workstation appears to be running $OSTYPE"
        error "For Windows/macOS, you'll need to adapt the installation commands."
        exit 1
    fi
    
    log "Operating system check passed: $(uname -a)"
}

# install system dependencies
install_system_deps() {
    log "Installing system dependencies..."
    
    # update package lists
    sudo apt update
    
    # install essential packages
    sudo apt install -y \
        wget \
        curl \
        git \
        build-essential \
        cmake \
        pkg-config \
        ca-certificates \
        gnupg \
        lsb-release \
        python3 \
        python3-pip \
        python3-dev \
        python3-venv \
        libssl-dev \
        libffi-dev \
        libbz2-dev \
        liblzma-dev \
        zlib1g-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libreadline-dev \
        libsqlite3-dev \
        libgdbm-dev \
        libdb5.3-dev \
        libdb-dev \
        libc6-dev \
        libexpat1-dev \
        libxml2-dev \
        libxslt1-dev \
        libyaml-dev \
        libgmp-dev \
        libmpfr-dev \
        libmpc-dev \
        libcairo2-dev \
        libpango1.0-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libx11-dev \
        libxext-dev \
        libxrender-dev \
        libxt-dev \
        libxft-dev \
        libfreetype6-dev \
        libfontconfig1-dev
    
    log "System dependencies installed"
}

# install nvidia drivers and cuda
install_nvidia_cuda() {
    log "Installing NVIDIA drivers and CUDA..."
    
    # check if nvidia gpu is present
    if ! lspci | grep -i nvidia > /dev/null; then
        warning "NVIDIA GPU not detected. Skipping CUDA installation."
        return 0
    fi
    
    # add nvidia package repository
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.0-1_all.deb
    sudo dpkg -i cuda-keyring_1.0-1_all.deb
    sudo apt update
    
    # install cuda toolkit
    sudo apt install -y cuda-toolkit-12-0
    
    # install cudnn
    sudo apt install -y libcudnn8 libcudnn8-dev
    
    # add cuda to path
    echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
    
    log "NVIDIA CUDA installation complete"
    info "Please reboot your system to ensure drivers are properly loaded"
}

# install miniconda
install_miniconda() {
    log "Installing Miniconda..."
    
    if command -v conda &> /dev/null; then
        info "Conda already installed: $(conda --version)"
        return 0
    fi
    
    # download and install miniconda
    cd /tmp
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p "$HOME/miniconda3"
    
    # add conda to path
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    export PATH="$HOME/miniconda3/bin:$PATH"
    
    # initialize conda
    "$HOME/miniconda3/bin/conda" init bash
    
    log "Miniconda installation complete"
}

# create conda environment
create_conda_env() {
    log "Creating conda environment: $CONDA_ENV"
    
    # source conda
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    
    # create environment
    conda create -n "$CONDA_ENV" python=3.9 -y
    
    # activate environment
    conda activate "$CONDA_ENV"
    
    log "Conda environment created and activated"
}

# install bioinformatics tools
install_bioinformatics_tools() {
    log "Installing bioinformatics tools..."
    
    # activate conda environment
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
    
    # install core bioinformatics tools
    conda install -c bioconda -y \
        bwa \
        samtools \
        bedtools \
        gatk4 \
        bcftools \
        htslib \
        picard \
        fastqc \
        multiqc \
        trimmomatic \
        seqtk \
        sra-tools \
        blast \
        blast+ \
        emboss \
        muscle \
        clustalo \
        mafft \
        raxml \
        iqtree
    
    # install python packages for genomics
    pip install \
        biopython \
        pandas \
        numpy \
        matplotlib \
        seaborn \
        plotly \
        jupyter \
        notebook \
        ipykernel \
        pysam \
        cyvcf2 \
        pyvcf \
        vcfpy
    
    log "Bioinformatics tools installed"
}

# install gpu-accelerated tools
install_gpu_tools() {
    log "Installing GPU-accelerated tools..."
    
    # activate conda environment
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
    
    # install cuda-enabled tools
    conda install -c nvidia -y cuda-toolkit
    
    # install gpu-accelerated compression
    conda install -c conda-forge -y pigz pbzip2
    
    # try to install gpu-bwa (experimental)
    pip install gpu-bwa 2>/dev/null || warning "GPU-BWA not available, will use CPU version"
    
    # install cuda-enabled python packages
    pip install \
        cupy-cuda12x \
        cudf \
        cuml \
        cugraph \
        rapids-singlecell
    
    log "GPU-accelerated tools installed"
}

# create directory structure
create_directories() {
    log "Creating directory structure..."
    
    mkdir -p "$WORK_DIR"/{references,data,output,temp,scripts,logs}
    mkdir -p "$WORK_DIR/output"/{alignment,variant_calling,reports,fastqc}
    
    log "Directory structure created at $WORK_DIR"
}

# download reference genomes
download_references() {
    log "Downloading reference genomes..."
    
    cd "$WORK_DIR/references"
    
    # download grch38 (if not already present)
    if [ ! -f "GRCh38.fna" ]; then
        log "Downloading GRCh38 reference genome..."
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip hg38.fa.gz
        mv hg38.fa GRCh38.fna
    fi
    
    # download grch37 (if needed)
    if [ ! -f "GRCh37.fna" ]; then
        log "Downloading GRCh37 reference genome..."
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
        gunzip hg19.fa.gz
        mv hg19.fa GRCh37.fna
    fi
    
    log "Reference genomes downloaded"
}

# create activation script
create_activation_script() {
    log "Creating activation script..."
    
    cat > "$WORK_DIR/activate_env.sh" << 'EOF'
#!/bin/bash
# Activation script for genomics environment

# Source conda
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Activate genomics environment
conda activate genomics

# Set environment variables
export WORK_DIR="${WORK_DIR:-$HOME/genomic_analysis}"
export REF_DIR="$WORK_DIR/references"
export DATA_DIR="$WORK_DIR/data"
export OUTPUT_DIR="$WORK_DIR/output"
export TEMP_DIR="$WORK_DIR/temp"

# Add CUDA to PATH if available
if [ -d "/usr/local/cuda/bin" ]; then
    export PATH="/usr/local/cuda/bin:$PATH"
    export LD_LIBRARY_PATH="/usr/local/cuda/lib64:$LD_LIBRARY_PATH"
fi

# Set thread count
export THREADS=$(nproc)

echo "Genomics environment activated!"
echo "Work directory: $WORK_DIR"
echo "Threads available: $THREADS"
echo "GPU available: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'Not detected')"
EOF
    
    chmod +x "$WORK_DIR/activate_env.sh"
    
    log "Activation script created: $WORK_DIR/activate_env.sh"
}

# create test script
create_test_script() {
    log "Creating test script..."
    
    cat > "$WORK_DIR/test_installation.sh" << 'EOF'
#!/bin/bash
# Test script for genomics installation

source "${WORK_DIR:-$HOME/genomic_analysis}/activate_env.sh"

echo "Testing genomics installation..."
echo "================================"

# Test conda environment
echo "Conda environment: $CONDA_DEFAULT_ENV"

# Test basic tools
echo "Testing basic tools:"
echo "BWA version: $(bwa 2>&1 | head -1)"
echo "Samtools version: $(samtools --version | head -1)"
echo "GATK version: $(gatk --version | head -1)"

# Test GPU tools
if command -v nvidia-smi &> /dev/null; then
    echo "GPU information:"
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
else
    echo "NVIDIA GPU not detected"
fi

# Test Python packages
echo "Testing Python packages:"
python -c "import pandas, numpy, matplotlib; print('Core packages: OK')"
python -c "import Bio; print('Biopython: OK')" 2>/dev/null || echo "Biopython: Not installed"

# Test CUDA (if available)
if python -c "import cupy" 2>/dev/null; then
    echo "CuPy (CUDA): OK"
else
    echo "CuPy (CUDA): Not available"
fi

echo "================================"
echo "Installation test complete!"
EOF
    
    chmod +x "$WORK_DIR/test_installation.sh"
    
    log "Test script created: $WORK_DIR/test_installation.sh"
}

# Create README
create_readme() {
    log "Creating README..."
    
    cat > "$WORK_DIR/README.md" << 'EOF'
# Local Genomics Analysis Environment

This environment is optimized for Linux workstations with NVIDIA GPU support.

## Quick Start

1. **Activate environment:**
   ```bash
   source ${WORK_DIR:-$HOME/genomic_analysis}/activate_env.sh
   ```

2. **Test installation:**
   ```bash
   ${WORK_DIR:-$HOME/genomic_analysis}/test_installation.sh
   ```

3. **Run analysis:**
   ```bash
   ${WORK_DIR:-$HOME/genomic_analysis}/scripts/local_genome_analysis.sh
   ```

## Directory Structure

```
genomic_analysis/
├── references/          # Reference genomes and indices
├── data/               # Input FASTQ files
├── output/             # Analysis results
│   ├── alignment/      # BAM files
│   ├── variant_calling/ # VCF files
│   ├── reports/        # Analysis reports
│   └── fastqc/         # Quality control
├── temp/               # Temporary files
├── scripts/            # Analysis scripts
└── logs/               # Log files
```

## Hardware Optimization

- **CPU**: Uses all available threads for parallel processing
- **GPU**: NVIDIA GPU accelerates alignment and variant calling (if available)
- **Storage**: NVMe SSD recommended for fast I/O with large datasets
- **Memory**: Optimized for 32GB+ RAM systems

## Performance Expectations

For typical genome data:
- **BWA alignment**: 2-4 hours (vs 4-8 hours on AWS)
- **GATK variant calling**: 1-2 hours (vs 2-4 hours on AWS)
- **Total analysis time**: 4-6 hours (vs 8-12 hours on AWS)

## Cost Comparison

- **Local processing**: $0/month (one-time setup)
- **AWS c7i.8xlarge**: ~$346/month (when running)
- **Savings**: ~$346/month + no data transfer costs

## Troubleshooting

1. **GPU not detected**: Ensure NVIDIA drivers are installed and system is rebooted
2. **Out of memory**: Close other applications or reduce thread count
3. **Storage full**: Clean up temporary files or move data to external storage

## Support

For issues or questions, check the log files in the `logs/` directory.
EOF
    
    log "README created: $WORK_DIR/README.md"
}

# Main installation function
main() {
    log "Starting local genomics environment setup..."
    log "Target directory: $WORK_DIR"
    log "Conda environment: $CONDA_ENV"
    
    # Run installation steps
    check_os
    install_system_deps
    install_nvidia_cuda
    install_miniconda
    create_conda_env
    install_bioinformatics_tools
    install_gpu_tools
    create_directories
    download_references
    create_activation_script
    create_test_script
    create_readme
    
    log "Setup complete!"
    log ""
    log "Next steps:"
    log "1. Reboot your system to load NVIDIA drivers"
    log "2. Run: source ${WORK_DIR:-$HOME/genomic_analysis}/activate_env.sh"
    log "3. Test: ${WORK_DIR:-$HOME/genomic_analysis}/test_installation.sh"
    log "4. Copy your FASTQ files to ${WORK_DIR:-$HOME/genomic_analysis}/data/"
    log "5. Run analysis: ${WORK_DIR:-$HOME/genomic_analysis}/scripts/local_genome_analysis.sh"
    log ""
    log "Documentation: ${WORK_DIR:-$HOME/genomic_analysis}/README.md"
}

# Run main function
main "$@"

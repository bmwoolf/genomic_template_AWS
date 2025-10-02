#!/bin/bash

# =============================================================================
# GENOME ANALYSIS WORKFLOW
# runs the actual genomic analysis pipeline
# requires: setup_bioconda.sh to be run first
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION - EDIT THESE FOR YOUR ENVIRONMENT
# =============================================================================
SAMPLE_NAME="${SAMPLE_NAME:-$(whoami)_genome}"  # defaults to your username
REFERENCE_GENOME="GRCh38"
THREADS="${THREADS:-$(nproc)}"  # defaults to all available CPU threads
GPU_ENABLED=true
WORK_DIR="${WORK_DIR:-$HOME/wgs_data}"  # defaults to ~/wgs_data
REF_DIR="$WORK_DIR/references"
DATA_DIR="$WORK_DIR/uncompressed"
OUTPUT_DIR="$WORK_DIR/output"
TEMP_DIR="$WORK_DIR/temp"

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

# check if environment is set up
check_environment() {
    log "Checking if environment is set up..."
    
    # check if conda is installed
    if ! command -v conda &> /dev/null; then
        error "Conda not found. Run './scripts/workstation/setup_bioconda.sh' first."
        exit 1
    fi
    
    # check if genomics environment exists
    if ! conda env list | grep -q "genomics"; then
        error "Genomics environment not found. Run './scripts/workstation/setup_bioconda.sh' first."
        exit 1
    fi
    
    # activate environment
    eval "$(conda shell.bash hook)"
    conda activate genomics
    
    # check if required tools are installed
    for tool in bwa samtools gatk bcftools multiqc; do
        if ! command -v $tool &> /dev/null; then
            error "$tool not found. Run './scripts/workstation/setup_bioconda.sh' first."
            exit 1
        fi
    done
    
    # check if reference genome exists
    if [ ! -f "$REF_DIR/GRCh38.fna" ]; then
        error "Reference genome not found. Run './scripts/workstation/setup_bioconda.sh' first."
        exit 1
    fi
    
    info "Environment check passed"
}

# check system requirements
check_requirements() {
    log "Checking system requirements..."
    
    # check available memory
    TOTAL_MEM=$(free -g | awk '/^Mem:/{print $2}')
    if [ $TOTAL_MEM -lt 32 ]; then
        warning "Less than 32GB RAM detected. You have ${TOTAL_MEM}GB. Consider closing other applications."
    else
        info "RAM check passed: ${TOTAL_MEM}GB available"
    fi
    
    # check available storage
    AVAILABLE_SPACE=$(df -BG "$WORK_DIR" | awk 'NR==2{print $4}' | sed 's/G//')
    if [ $AVAILABLE_SPACE -lt 1200 ]; then
        error "Less than 1.2TB free space available. Need more storage for 485GB genome analysis."
        exit 1
    else
        info "Storage check passed: ${AVAILABLE_SPACE}GB available"
    fi
    
    # check GPU
    if command -v nvidia-smi &> /dev/null; then
        GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
        GPU_MEM_GB=$((GPU_MEM / 1024))
        info "GPU detected: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
        info "GPU Memory: ${GPU_MEM_GB}GB"
        
        if [ $GPU_MEM_GB -lt 8 ]; then
            warning "GPU has less than 8GB VRAM. Some GPU acceleration may be limited."
        else
            info "GPU VRAM check passed: ${GPU_MEM_GB}GB available"
        fi
    else
        warning "NVIDIA GPU not detected. GPU acceleration disabled."
        GPU_ENABLED=false
    fi
    
    # check CPU
    CPU_CORES=$(nproc)
    info "CPU cores detected: $CPU_CORES"
    
    log "System requirements check complete"
}

# process FASTQ files
process_fastq() {
    log "Processing FASTQ files..."
    
    mkdir -p "$DATA_DIR" "$OUTPUT_DIR" "$TEMP_DIR"
    cd "$DATA_DIR"
    
    # check if FASTQ files exist
    if [ ! -f *.gz ] && [ ! -f *.fastq ]; then
        error "No FASTQ files found in $DATA_DIR"
        error "Please place your FASTQ files (.gz or .fastq) in $DATA_DIR"
        exit 1
    fi
    
    # decompress FASTQ files if needed
    if ls *.gz 1> /dev/null 2>&1; then
        log "Decompressing FASTQ files..."
        for file in *.gz; do
            if [ -f "$file" ]; then
                log "Decompressing $file..."
                pigz -d -p $THREADS "$file"
            fi
        done
    fi
    
    # quality control
    log "Running quality control..."
    mkdir -p "$OUTPUT_DIR/fastqc"
    fastqc *.fastq -o "$OUTPUT_DIR/fastqc" -t $THREADS
    
    # generate multiQC report
    log "Generating MultiQC report..."
    multiqc "$OUTPUT_DIR/fastqc" -o "$OUTPUT_DIR" --force
    
    log "FASTQ processing complete"
}

# BWA alignment
run_alignment() {
    log "Running BWA alignment..."
    
    cd "$DATA_DIR"
    
    # find FASTQ files
    R1_FILES=($(ls *R1*.fastq 2>/dev/null || ls *_1.fastq 2>/dev/null || echo ""))
    R2_FILES=($(ls *R2*.fastq 2>/dev/null || ls *_2.fastq 2>/dev/null || echo ""))
    
    if [ ${#R1_FILES[@]} -eq 0 ]; then
        error "No R1 FASTQ files found. Expected files matching *R1*.fastq or *_1.fastq"
        exit 1
    fi
    
    # create output directory
    mkdir -p "$OUTPUT_DIR/alignment"
    
    # run alignment for each sample
    for i in "${!R1_FILES[@]}"; do
        R1_FILE="${R1_FILES[$i]}"
        R2_FILE="${R2_FILES[$i]}"
        
        log "Aligning $R1_FILE and $R2_FILE..."
        
        # extract lane info from filename (e.g., L001 from filename)
        LANE=$(echo "$R1_FILE" | grep -oP 'L\d+' || echo "L00${i}")
        READ_GROUP="@RG\tID:${LANE}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:lib1\tPU:${LANE}"
        
        bwa mem -t $THREADS -R "$READ_GROUP" \
            "$REF_DIR/GRCh38.fna" \
            "$R1_FILE" "$R2_FILE" \
            | samtools view -bS - \
            | samtools sort -@ $THREADS -o "$OUTPUT_DIR/alignment/${SAMPLE_NAME}_${i}.bam" -
        
        # index BAM file
        samtools index "$OUTPUT_DIR/alignment/${SAMPLE_NAME}_${i}.bam"
        
        log "Alignment complete for sample $i"
    done
    
    # merge BAM files if multiple samples
    if [ ${#R1_FILES[@]} -gt 1 ]; then
        log "Merging BAM files..."
        samtools merge -@ $THREADS \
            "$OUTPUT_DIR/alignment/${SAMPLE_NAME}_merged.bam" \
            "$OUTPUT_DIR/alignment/${SAMPLE_NAME}_"*.bam
        
        samtools index "$OUTPUT_DIR/alignment/${SAMPLE_NAME}_merged.bam"
    fi
    
    log "Alignment complete"
}

# GATK variant calling
run_variant_calling() {
    log "Running GATK variant calling..."
    
    cd "$OUTPUT_DIR"
    mkdir -p "variant_calling"
    
    # determine input BAM file
    if [ -f "alignment/${SAMPLE_NAME}_merged.bam" ]; then
        INPUT_BAM="alignment/${SAMPLE_NAME}_merged.bam"
    else
        INPUT_BAM="alignment/${SAMPLE_NAME}_0.bam"
    fi
    
    # mark duplicates
    log "Marking duplicates..."
    gatk MarkDuplicates \
        -I "$INPUT_BAM" \
        -O "alignment/${SAMPLE_NAME}_dedup.bam" \
        -M "alignment/${SAMPLE_NAME}_metrics.txt" \
        --TMP_DIR "$TEMP_DIR"
    
    samtools index "alignment/${SAMPLE_NAME}_dedup.bam"
    
    # base quality score recalibration
    log "Running BQSR..."
    gatk BaseRecalibrator \
        -R "$REF_DIR/GRCh38.fna" \
        -I "alignment/${SAMPLE_NAME}_dedup.bam" \
        --known-sites "$REF_DIR/gatk_resources/Homo_sapiens_assembly38.known_indels.vcf.gz" \
        --known-sites "$REF_DIR/gatk_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
        --known-sites "$REF_DIR/gatk_resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \
        -O "variant_calling/${SAMPLE_NAME}_recal.table" \
        --tmp-dir "$TEMP_DIR"
    
    # apply BQSR
    log "Applying BQSR..."
    gatk ApplyBQSR \
        -R "$REF_DIR/GRCh38.fna" \
        -I "alignment/${SAMPLE_NAME}_dedup.bam" \
        -bqsr "variant_calling/${SAMPLE_NAME}_recal.table" \
        -O "alignment/${SAMPLE_NAME}_recal.bam" \
        --tmp-dir "$TEMP_DIR"
    
    samtools index "alignment/${SAMPLE_NAME}_recal.bam"
    
    # variant calling
    log "Running HaplotypeCaller..."
    gatk HaplotypeCaller \
        -R "$REF_DIR/GRCh38.fna" \
        -I "alignment/${SAMPLE_NAME}_recal.bam" \
        -O "variant_calling/${SAMPLE_NAME}_raw.vcf.gz" \
        --native-pair-hmm-threads $THREADS \
        --tmp-dir "$TEMP_DIR"
    
    # variant filtering
    log "Filtering variants..."
    gatk VariantFiltration \
        -R "$REF_DIR/GRCh38.fna" \
        -V "variant_calling/${SAMPLE_NAME}_raw.vcf.gz" \
        -O "variant_calling/${SAMPLE_NAME}_filtered.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \
        --filter-expression "FS > 60.0" --filter-name "FS60" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
    
    log "Variant calling complete"
}

# generate final reports
generate_reports() {
    log "Generating final reports..."
    
    cd "$OUTPUT_DIR"
    mkdir -p "reports"
    
    # variant statistics
    log "Generating variant statistics..."
    bcftools stats "variant_calling/${SAMPLE_NAME}_filtered.vcf.gz" > "reports/${SAMPLE_NAME}_variant_stats.txt"
    
    # coverage analysis
    log "Generating coverage analysis..."
    gatk CollectWgsMetrics \
        -R "$REF_DIR/GRCh38.fna" \
        -I "alignment/${SAMPLE_NAME}_recal.bam" \
        -O "reports/${SAMPLE_NAME}_coverage.txt"
    
    # create summary report
    log "Creating summary report..."
    cat > "reports/${SAMPLE_NAME}_summary.txt" << EOF
Genomic Analysis Summary
========================
Sample: $SAMPLE_NAME
Reference: $REFERENCE_GENOME
Analysis Date: $(date)
Workstation: Linux with NVIDIA GPU support

Input Data:
- FASTQ files: $(ls "$DATA_DIR"/*.fastq 2>/dev/null | wc -l) files
- Total size: $(du -sh "$DATA_DIR" 2>/dev/null | cut -f1)

Output Files:
- BAM file: alignment/${SAMPLE_NAME}_recal.bam
- VCF file: variant_calling/${SAMPLE_NAME}_filtered.vcf.gz
- Coverage report: reports/${SAMPLE_NAME}_coverage.txt
- Variant stats: reports/${SAMPLE_NAME}_variant_stats.txt

Performance Notes:
- CPU threads used: $THREADS
- GPU acceleration: $GPU_ENABLED
- Analysis completed: $(date)
EOF
    
    log "Reports generated successfully"
}

# cleanup temporary files
cleanup() {
    log "Cleaning up temporary files..."
    
    # remove temporary directory contents
    rm -rf "$TEMP_DIR"/*
    
    # compress large intermediate files
    log "Compressing intermediate files..."
    pigz -p $THREADS "$OUTPUT_DIR/alignment/${SAMPLE_NAME}_dedup.bam" 2>/dev/null || true
    
    log "Cleanup complete"
}

# main execution
main() {
    log "==================================================================="
    log "Starting genomic analysis workflow"
    log "==================================================================="
    log "Sample: $SAMPLE_NAME"
    log "Reference: $REFERENCE_GENOME"
    log "Threads: $THREADS"
    log "GPU enabled: $GPU_ENABLED"
    
    # create directory structure
    mkdir -p "$WORK_DIR" "$DATA_DIR" "$OUTPUT_DIR" "$TEMP_DIR" "$OUTPUT_DIR/reports"
    
    # run analysis steps
    check_environment
    check_requirements
    process_fastq
    run_alignment
    run_variant_calling
    generate_reports
    cleanup
    
    log "==================================================================="
    log "CPU Variant Calling Complete! Results available in $OUTPUT_DIR"
    log "==================================================================="
    log "Summary report: $OUTPUT_DIR/reports/${SAMPLE_NAME}_summary.txt"
    log "Final VCF: $OUTPUT_DIR/variant_calling/${SAMPLE_NAME}_filtered.vcf.gz"
    log ""
    log "Next: Run annotation and prioritization:"
    log "  ./scripts/workstation/4_annotate_and_prioritize.sh"
}

# run main function
main "$@"

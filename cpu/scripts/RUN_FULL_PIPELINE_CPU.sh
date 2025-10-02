#!/bin/bash

# =============================================================================
# MASTER GENOMIC ANALYSIS PIPELINE
# runs complete workflow: Alignment → Variant Calling → Annotation → Prioritization
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# CONFIGURATION - EDIT THESE FOR YOUR ENVIRONMENT
# =============================================================================
export SAMPLE_NAME="${SAMPLE_NAME:-$(whoami)_genome}"
export WORK_DIR="${WORK_DIR:-$HOME/wgs_data}"
export THREADS="${THREADS:-$(nproc)}"

# colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log "==================================================================="
log "COMPLETE CPU-BASED GENOMIC ANALYSIS PIPELINE"
log "==================================================================="
log ""

# step 1: set up Linux-based computer environment
log "step 1: set up Linux-based computer environment"
info "this includes: installing system dependencies, conda, and bioinformatics tools"
"$SCRIPT_DIR/1_setup_computer.sh"

# step 2: set up bioconda environment
log "step 2: set up bioconda environment for running end to end genomic analysis pipeline"
info "this includes: installing bioinformatics tools via conda"
"$SCRIPT_DIR/2_setup_bioconda.sh"

# step 3: run FASTQ -> VCF
log "step 3: running complete genomic analysis pipeline (FASTQ -> VCF)"
info "this includes: QC → Alignment → BQSR → Variant Calling → Filtering"
"$SCRIPT_DIR/3_analyze_genome_cpu.sh"

# step 4: run VCF -> analysis via annotation and prioritization
log "step 4: running analysis via annotation and prioritization"
"$SCRIPT_DIR/4_annotate_and_prioritize.sh"

log ""
log "==================================================================="
log "COMPLETE PIPELINE FINISHED!"
log "==================================================================="
log ""
log "Sample: $SAMPLE_NAME"
log "Output directory: $WORK_DIR/output"
log ""
log "Key output files:"
log "  - VCF: $WORK_DIR/output/variant_calling/${SAMPLE_NAME}_filtered.vcf.gz"
log "  - Annotated VCF: $WORK_DIR/output/annotated/${SAMPLE_NAME}_vep_annotated.vcf.gz"
log "  - Prioritized variants: $WORK_DIR/output/annotated/${SAMPLE_NAME}_prioritized_variants.tsv"
log "  - Structural variants: $WORK_DIR/output/structural_variants/manta/"
log "  - CNVs: $WORK_DIR/output/structural_variants/cnvnator/"
log "" 
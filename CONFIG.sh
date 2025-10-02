#!/bin/bash

# =============================================================================
# GENOMIC ANALYSIS PIPELINE - CENTRAL CONFIGURATION
# edit these settings for your environment
# =============================================================================

# sample name (will be used in output filenames)
# default: uses your system username + "_genome"
# example: "john_doe_genome" or "patient_001"
export SAMPLE_NAME="${SAMPLE_NAME:-$(whoami)_genome}"

# working directory (where all data and results will be stored)
# default: ~/wgs_data (in your home directory)
# requirement: needs ~1.5TB+ free space for whole genome analysis
export WORK_DIR="${WORK_DIR:-$HOME/wgs_data}"

# number of CPU threads to use
# default: uses all available CPU threads
# you can manually set this to a lower number if needed
export THREADS=$(nproc)

# reference genome version
export REFERENCE_GENOME="GRCh38"

# GPU acceleration (set to false if you don't have a GPU)
export GPU_ENABLED=true

# =============================================================================
# ADVANCED SETTINGS (usually don't need to change)
# =============================================================================

# directory paths (relative to WORK_DIR)
export REF_DIR="$WORK_DIR/references"
export DATA_DIR="$WORK_DIR/uncompressed"
export OUTPUT_DIR="$WORK_DIR/output"
export TEMP_DIR="$WORK_DIR/temp"

# conda environment name
export CONDA_ENV="genomics"

# =============================================================================
# USAGE
# =============================================================================
# to use this configuration, either:
# 1. edit the values above and source this file before running scripts:
#    source CONFIG.sh
#    ./scripts/workstation/3_analyze_genome_cpu.sh
#
# 2. or set environment variables directly:
#    export SAMPLE_NAME="my_sample"
#    export WORK_DIR="/mnt/bigdrive/genomics"
#    ./scripts/workstation/RUN_FULL_PIPELINE.sh
# =============================================================================

echo "Configuration loaded:"
echo "  Sample: $SAMPLE_NAME"
echo "  Work directory: $WORK_DIR"
echo "  Threads: $THREADS"
echo "  Reference: $REFERENCE_GENOME" 
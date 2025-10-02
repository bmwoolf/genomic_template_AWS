#!/bin/bash

# =============================================================================
# POST-VCF ANNOTATION AND PRIORITIZATION
# runs SnpEff annotation, variant prioritization, and SV/CNV detection
# prerequisites: VCF file from step 3
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION - SOURCE CONFIG.SH OR USE DEFAULTS
# =============================================================================
# unset any previously set variables to force fresh config
unset SAMPLE_NAME WORK_DIR

# try to source CONFIG.sh from script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "$SCRIPT_DIR/../../CONFIG.sh" ]; then
    source "$SCRIPT_DIR/../../CONFIG.sh"
elif [ -f "$SCRIPT_DIR/../CONFIG.sh" ]; then
    source "$SCRIPT_DIR/../CONFIG.sh"
elif [ -f "CONFIG.sh" ]; then
    source CONFIG.sh
fi

# set defaults if not already set by CONFIG.sh
SAMPLE_NAME="${SAMPLE_NAME:-$(whoami)_genome}"
REFERENCE_GENOME="GRCh38"
THREADS="${THREADS:-$(nproc)}"
WORK_DIR="${WORK_DIR:-$HOME/wgs_data}"
REF_DIR="$WORK_DIR/references"
OUTPUT_DIR="$WORK_DIR/output"
VCF_DIR="$OUTPUT_DIR/variant_calling"
ANNOTATED_DIR="$OUTPUT_DIR/annotated"
SV_DIR="$OUTPUT_DIR/structural_variants"

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

# create annotation environment if needed
setup_annotation_environment() {
    log "Setting up annotation environment..."
    
    eval "$(conda shell.bash hook)"
    
    # check if annotation environment exists
    if ! conda env list | grep -q "^annotation "; then
        log "Creating new conda environment for annotation tools..."
        conda create -n annotation python=3.10 -y
    fi
    
    conda activate annotation
    
    # check if tools are installed
    if ! command -v bcftools &> /dev/null; then
        log "Installing core annotation tools (5-10 minutes)..."
        conda install -c conda-forge -c bioconda \
            bcftools \
            samtools \
            tabix \
            htslib \
            delly \
            cnvkit \
            -y
        
        info "Core annotation tools installed successfully"
    else
        info "Annotation tools already installed"
    fi
}

# install SnpEff for variant annotation (reliable Java-based alternative to VEP)
install_snpeff() {
    log "Installing SnpEff for variant annotation..."
    
    eval "$(conda shell.bash hook)"
    conda activate annotation
    
    # check if snpEff is already installed
    if command -v snpEff &> /dev/null; then
        info "SnpEff already installed"
        return 0
    fi
    
    log "Installing SnpEff via conda (2-3 minutes)..."
    conda install -c bioconda snpeff -y || {
        warning "SnpEff installation failed"
        return 1
    }
    
    # download GRCh38 database
    log "Downloading GRCh38 annotation database (this may take 5-10 minutes)..."
    snpEff download -v GRCh38.mane.1.0.refseq || {
        warning "Using fallback database..."
        snpEff download -v GRCh38.86 || {
            warning "Database download failed - will download on first use"
            return 1
        }
    }
    
    info "SnpEff installed successfully"
    return 0
}

# run SnpEff annotation
run_snpeff_annotation() {
    log "Running SnpEff annotation..."
    
    eval "$(conda shell.bash hook)"
    conda activate annotation
    
    mkdir -p "$ANNOTATED_DIR"
    
    local INPUT_VCF="$VCF_DIR/${SAMPLE_NAME}_filtered.vcf.gz"
    local OUTPUT_VCF="$ANNOTATED_DIR/${SAMPLE_NAME}_annotated.vcf"
    local OUTPUT_HTML="$ANNOTATED_DIR/${SAMPLE_NAME}_snpeff_summary.html"
    local OUTPUT_TSV="$ANNOTATED_DIR/${SAMPLE_NAME}_annotated.tsv"
    
    if [ ! -f "$INPUT_VCF" ]; then
        error "Input VCF not found: $INPUT_VCF"
        exit 1
    fi
    
    if command -v snpEff &> /dev/null; then
        log "Annotating variants with SnpEff..."
        
        # try MANE database first (most current), fall back to 86
        local DATABASE="GRCh38.mane.1.0.refseq"
        if ! snpEff databases | grep -q "$DATABASE"; then
            DATABASE="GRCh38.86"
        fi
        
        log "Using database: $DATABASE"
        
        # set Java memory options (use up to 16GB for faster processing)
        export _JAVA_OPTIONS="-Xms8g -Xmx16g"
        
        # run SnpEff annotation
        snpEff ann \
            -v \
            -stats "$OUTPUT_HTML" \
            -csvStats "${OUTPUT_HTML%.html}.csv" \
            "$DATABASE" \
            "$INPUT_VCF" > "$OUTPUT_VCF" || {
                warning "SnpEff annotation failed, falling back to basic annotation"
                run_basic_annotation
                return 1
            }
        
        # unset to avoid affecting other tools
        unset _JAVA_OPTIONS
        
        # compress and index VCF
        bgzip -f "$OUTPUT_VCF"
        tabix -p vcf "$OUTPUT_VCF.gz"
        
        # extract to TSV for easy viewing
        log "Creating TSV output..."
        echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENE\tEFFECT\tIMPACT\tANN" > "$OUTPUT_TSV"
        
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/ANN\n' \
            "$OUTPUT_VCF.gz" | \
        awk -F'\t' '{
            # parse ANN field (pipe-delimited)
            split($8, ann, "|");
            gene = ann[4];
            effect = ann[2];
            impact = ann[3];
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
                $1, $2, $3, $4, $5, $6, $7, gene, effect, impact, $8;
        }' >> "$OUTPUT_TSV"
        
        log "SnpEff annotation complete"
        info "Summary report: $OUTPUT_HTML"
        info "Annotated VCF: $OUTPUT_VCF.gz"
        info "TSV export: $OUTPUT_TSV"
    else
        warning "SnpEff not available, using basic annotation"
        run_basic_annotation
    fi
}

# basic annotation fallback
run_basic_annotation() {
    log "Running basic variant annotation..."
    
    local INPUT_VCF="$VCF_DIR/${SAMPLE_NAME}_filtered.vcf.gz"
    local OUTPUT_TSV="$ANNOTATED_DIR/${SAMPLE_NAME}_basic_annotation.tsv"
    
    # extract variant information to TSV
    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$OUTPUT_TSV"
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' \
        "$INPUT_VCF" >> "$OUTPUT_TSV"
    
    log "Basic annotation complete: $OUTPUT_TSV"
}

# run Delly for structural variant detection (replaces Manta)
run_delly_sv() {
    log "Running Delly for structural variant detection..."
    
    eval "$(conda shell.bash hook)"
    conda activate annotation
    
    mkdir -p "$SV_DIR/delly"
    
    local BAM_FILE="$OUTPUT_DIR/alignment/${SAMPLE_NAME}_recal.bam"
    
    if [ ! -f "$BAM_FILE" ]; then
        # try dedup_rg.bam instead
        BAM_FILE="$OUTPUT_DIR/alignment/${SAMPLE_NAME}_dedup_rg.bam"
        if [ ! -f "$BAM_FILE" ]; then
            error "BAM file not found"
            return 1
        fi
    fi
    
    # check if BAM is indexed
    if [ ! -f "$BAM_FILE.bai" ]; then
        log "Indexing BAM file..."
        samtools index "$BAM_FILE"
    fi
    
    # run Delly for different SV types IN PARALLEL (4x speedup!)
    log "Detecting all SV types in parallel (DEL, DUP, INV, BND)..."
    
    # run all 4 SV types in background
    delly call -t DEL -g "$REF_DIR/GRCh38.fna" -o "$SV_DIR/delly/${SAMPLE_NAME}.del.bcf" "$BAM_FILE" &
    local PID_DEL=$!
    
    delly call -t DUP -g "$REF_DIR/GRCh38.fna" -o "$SV_DIR/delly/${SAMPLE_NAME}.dup.bcf" "$BAM_FILE" &
    local PID_DUP=$!
    
    delly call -t INV -g "$REF_DIR/GRCh38.fna" -o "$SV_DIR/delly/${SAMPLE_NAME}.inv.bcf" "$BAM_FILE" &
    local PID_INV=$!
    
    delly call -t BND -g "$REF_DIR/GRCh38.fna" -o "$SV_DIR/delly/${SAMPLE_NAME}.bnd.bcf" "$BAM_FILE" &
    local PID_BND=$!
    
    # wait for all to complete
    log "Waiting for DEL detection (PID: $PID_DEL)..."
    wait $PID_DEL && info "DEL complete" || warning "DEL failed"
    
    log "Waiting for DUP detection (PID: $PID_DUP)..."
    wait $PID_DUP && info "DUP complete" || warning "DUP failed"
    
    log "Waiting for INV detection (PID: $PID_INV)..."
    wait $PID_INV && info "INV complete" || warning "INV failed"
    
    log "Waiting for BND detection (PID: $PID_BND)..."
    wait $PID_BND && info "BND complete" || warning "BND failed"
    
    log "All Delly SV detection complete"
    
    # merge all SV types
    log "Merging SV calls..."
    bcftools concat \
        "$SV_DIR/delly/${SAMPLE_NAME}.del.bcf" \
        "$SV_DIR/delly/${SAMPLE_NAME}.dup.bcf" \
        "$SV_DIR/delly/${SAMPLE_NAME}.inv.bcf" \
        "$SV_DIR/delly/${SAMPLE_NAME}.bnd.bcf" \
        -a -O v -o "$SV_DIR/delly/${SAMPLE_NAME}.sv.vcf"
    
    # compress and index
    bgzip -f "$SV_DIR/delly/${SAMPLE_NAME}.sv.vcf"
    tabix -p vcf "$SV_DIR/delly/${SAMPLE_NAME}.sv.vcf.gz"
    
    log "Delly SV detection complete"
    log "Results in: $SV_DIR/delly/${SAMPLE_NAME}.sv.vcf.gz"
}

# run CNVkit for copy number variant detection (replaces CNVnator)
run_cnvkit() {
    log "Running CNVkit for CNV detection..."
    
    eval "$(conda shell.bash hook)"
    conda activate annotation
    
    mkdir -p "$SV_DIR/cnvkit"
    
    local BAM_FILE="$OUTPUT_DIR/alignment/${SAMPLE_NAME}_recal.bam"
    
    if [ ! -f "$BAM_FILE" ]; then
        BAM_FILE="$OUTPUT_DIR/alignment/${SAMPLE_NAME}_dedup_rg.bam"
        if [ ! -f "$BAM_FILE" ]; then
            error "BAM file not found"
            return 1
        fi
    fi
    
    # CNVkit batch command (creates reference and calls CNVs)
    log "Running CNVkit batch analysis..."
    cnvkit.py batch \
        "$BAM_FILE" \
        --fasta "$REF_DIR/GRCh38.fna" \
        --output-dir "$SV_DIR/cnvkit" \
        --diagram \
        --scatter \
        -p $THREADS
    
    log "CNVkit analysis complete"
    log "Results in: $SV_DIR/cnvkit/"
}

# create VPOT-style prioritization (simplified, no external tool needed)
run_variant_prioritization() {
    log "Running variant prioritization..."
    
    eval "$(conda shell.bash hook)"
    conda activate annotation
    
    local INPUT_VCF="$VCF_DIR/${SAMPLE_NAME}_filtered.vcf.gz"
    local OUTPUT_TSV="$ANNOTATED_DIR/${SAMPLE_NAME}_prioritized_variants.tsv"
    
    log "Extracting high-quality variants..."
    
    # extract high-quality variants with detailed information
    cat > "$OUTPUT_TSV" << 'EOF'
# Prioritized Variant Report
# Filtered for: QUAL >= 30, FILTER = PASS
# Columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
EOF
    
    bcftools view -f PASS "$INPUT_VCF" | \
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' | \
        awk '$6 >= 30' >> "$OUTPUT_TSV"
    
    log "Prioritization complete: $OUTPUT_TSV"
    
    # count variants
    local TOTAL_VARIANTS=$(zcat "$INPUT_VCF" | grep -v "^#" | wc -l)
    local HIGH_QUALITY=$(grep -v "^#" "$OUTPUT_TSV" | wc -l)
    
    info "Total variants: $TOTAL_VARIANTS"
    info "High-quality variants (PASS, QUAL>=30): $HIGH_QUALITY"
}

# generate summary report
generate_summary() {
    log "Generating annotation summary..."
    
    local VARIANT_COUNT=$(zcat "$VCF_DIR/${SAMPLE_NAME}_filtered.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l || echo "Unknown")
    local HIGH_QUAL_COUNT=$(grep -v "^#" "$ANNOTATED_DIR/${SAMPLE_NAME}_prioritized_variants.tsv" 2>/dev/null | wc -l || echo "Unknown")
    
    cat > "$ANNOTATED_DIR/annotation_summary.txt" << EOF
Annotation and Prioritization Summary
=====================================
Sample: $SAMPLE_NAME
Analysis Date: $(date)

Input Files:
- VCF: $VCF_DIR/${SAMPLE_NAME}_filtered.vcf.gz
- Total variants: $VARIANT_COUNT
- High-quality variants (PASS, QUAL>=30): $HIGH_QUAL_COUNT

Analysis Tools:
- SnpEff: Genomic variant annotation (genes, effects, impacts)
- Delly: Structural variant detection (DEL, DUP, INV, BND)
- CNVkit: Copy number variant detection

Output Files:
- Annotated VCF: $ANNOTATED_DIR/${SAMPLE_NAME}_annotated.vcf.gz
- Annotation TSV: $ANNOTATED_DIR/${SAMPLE_NAME}_annotated.tsv
- Annotation report: $ANNOTATED_DIR/${SAMPLE_NAME}_snpeff_summary.html
- Prioritized variants: $ANNOTATED_DIR/${SAMPLE_NAME}_prioritized_variants.tsv
- Delly SVs: $SV_DIR/delly/${SAMPLE_NAME}.sv.vcf.gz
- CNVkit results: $SV_DIR/cnvkit/

Next Steps:
1. Open ${SAMPLE_NAME}_snpeff_summary.html in browser for annotation overview
2. Review prioritized variants (QUAL >= 30, FILTER = PASS)
3. Check annotated TSV for gene names, effects, and impact (HIGH/MODERATE/LOW)
4. Examine structural variants from Delly
5. Review CNV calls from CNVkit
6. Upload to VarSome.com or Franklin for additional clinical interpretation

Key Metrics:
- Coverage: Check reports/${SAMPLE_NAME}_coverage.txt
- Focus on HIGH and MODERATE impact variants
- Look for loss-of-function variants (stop_gained, frameshift_variant)
- Check missense variants in known disease genes

Notes:
- SnpEff provides: gene names, transcript IDs, variant effects, impact predictions
- IMPACT levels: HIGH (protein-truncating), MODERATE (missense), LOW (synonymous), MODIFIER (intergenic)
- Delly detects SVs >50bp (deletions, duplications, inversions, translocations)
- CNVkit provides CNV calls with visualization plots
- For clinical use, consult with a genetic counselor
EOF
    
    log "Summary saved to: $ANNOTATED_DIR/annotation_summary.txt"
}

# main execution
main() {
    log "==================================================================="
    log "Starting post-VCF annotation and prioritization"
    log "==================================================================="
    
    # setup environment
    setup_annotation_environment
    
    # install and run SnpEff annotation
    if install_snpeff; then
        run_snpeff_annotation || warning "SnpEff annotation failed - basic annotation used"
    else
        run_basic_annotation
    fi
    
    # run variant prioritization
    run_variant_prioritization
    
    # run SV detection with Delly
    run_delly_sv || warning "Delly SV detection failed - continuing anyway"
    
    # run CNV detection with CNVkit
    run_cnvkit || warning "CNVkit analysis failed - continuing anyway"
    
    # generate summary
    generate_summary
    
    log "==================================================================="
    log "Annotation and prioritization complete!"
    log "==================================================================="
    log "Summary: $ANNOTATED_DIR/annotation_summary.txt"
    log "Prioritized variants: $ANNOTATED_DIR/${SAMPLE_NAME}_prioritized_variants.tsv"
    log "Annotated TSV: $ANNOTATED_DIR/${SAMPLE_NAME}_annotated.tsv"
    log "Annotation report: $ANNOTATED_DIR/${SAMPLE_NAME}_snpeff_summary.html"
    log "Structural variants: $SV_DIR/delly/${SAMPLE_NAME}.sv.vcf.gz"
    log "CNV results: $SV_DIR/cnvkit/"
    log ""
    log "🎉 Your genome analysis is complete!"
}

# run main function
main "$@" 
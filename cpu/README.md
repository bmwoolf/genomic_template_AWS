# Local Workstation Analysis Scripts (Linux)

Complete end-to-end CPU-based genomic analysis pipeline for running on your own Linux workstation.

## Pipeline Overview

```
input: FASTQ files
    ↓
1. quality control (FastQC/MultiQC)
    ↓
2. alignment (BWA-MEM) → BAM
    ↓
3. mark duplicates (Picard)
    ↓
4. base quality score recalibration (GATK BQSR)
    ↓
5. variant calling (HaplotypeCaller)
    ↓
6. output VCF!
    ↓
7. variant filtering (GATK VariantFiltration)
    ↓
8. variant consequences (SnpEff)
    ↓
9. variant prioritization (bcftools + awk, QUAL>=30)
    ↓
10. structural variant detection (Delly: DEL/DUP/INV/BND)
    ↓
11. copy number variant detection (CNVkit)
    ↓
output: lots of information about your genome!
```

## Quick Start

### Step 1: Configure for Your Environment
Edit `CONFIG.sh`:
```bash
export SAMPLE_NAME="your_sample_name"  # or leave default: $(whoami)_genome
export WORK_DIR="$HOME/wgs_data"       # or your preferred directory
export THREADS=$(nproc)                # or set manually, e.g., 24
export REFERENCE_GENOME="GRCh38"
```

### Step 2: Run Complete Pipeline
```bash
source CONFIG.sh  # load your configuration
./scripts/RUN_FULL_PIPELINE_CPU.sh
```

## Output Files

### Variant Calling Outputs
```
$WORK_DIR/output/  (default: ~/wgs_data/output/)
├── alignment/
│   ├── {SAMPLE_NAME}_merged.bam          # Merged alignment
│   ├── {SAMPLE_NAME}_dedup_rg.bam        # Deduplicated with read groups
│   └── {SAMPLE_NAME}_recal.bam           # BQSR-recalibrated (final BAM)
│
├── variant_calling/
│   ├── {SAMPLE_NAME}_raw.vcf.gz          # Raw variants
│   └── {SAMPLE_NAME}_filtered.vcf.gz     # Filtered variants (FINAL VCF)
│
├── reports/
|   ├── {SAMPLE_NAME}_summary.txt          # Analysis summary
|   ├── {SAMPLE_NAME}_variant_stats.txt    # Variant statistics
|   └── {SAMPLE_NAME}_coverage.txt         # Coverage metrics
|
├── annotated/
│   ├── {SAMPLE_NAME}_annotated.vcf.gz           # SnpEff annotated VCF
│   ├── {SAMPLE_NAME}_annotated.tsv              # Gene/effect/impact table
│   ├── {SAMPLE_NAME}_snpeff_summary.html        # Interactive HTML report
│   ├── {SAMPLE_NAME}_prioritized_variants.tsv   # High-quality variants (TOP FILE)
│   └── annotation_summary.txt                    # Summary report
│
└── structural_variants/
    ├── delly/
    │   └── {SAMPLE_NAME}.sv.vcf.gz              # Structural variants (DEL/DUP/INV/BND)
    └── cnvkit/
        └── *.cns, *.cnr, *.pdf                  # Copy number variants + plots
```

Note: `{SAMPLE_NAME}` will be replaced with your configured sample name.

## Time Estimates

| Step | Estimated Time |
|------|---------------|
| Quality Control | 30 min |
| Alignment (BWA) | 12 hours |
| Mark Duplicates | 1 hour |
| BQSR | 2-3 hours |
| Variant Calling | 3-4 hours |
| Filtering | 30 min |
| SnpEff Annotation | 30-45 min |
| Variant Prioritization | 1 min |
| Delly SV (4 types) | 1-2 hours |
| CNVkit CNV | 30 min |
| **TOTAL** | **~21-25 hours** |

## Interpreting Results

### Result files
1. **`{SAMPLE_NAME}_snpeff_summary.html`**: START HERE (open in browser)
   - Interactive report with variant summary statistics
   - Charts showing variant types, effects, and impacts
   - Gene-level and chromosome-level distributions

2. **`{SAMPLE_NAME}_prioritized_variants.tsv`**: HIGH-QUALITY VARIANTS
   - 4.5M+ filtered variants (PASS, QUAL≥30)
   - Use for downstream analysis

3. **`{SAMPLE_NAME}_annotated.tsv`**: ANNOTATED VARIANTS
   - Gene names, variant effects (missense, stop_gained, etc.)
   - Impact levels: HIGH (protein-truncating), MODERATE (missense), LOW (synonymous)
   - Focus on HIGH and MODERATE impact variants

4. **`delly/{SAMPLE_NAME}.sv.vcf.gz`**: STRUCTURAL VARIANTS
   - Deletions (DEL), duplications (DUP), inversions (INV), translocations (BND)
   - Check for known disease-associated SVs >50bp

5. **`cnvkit/`**: COPY NUMBER VARIANTS
   - Gene duplications/deletions with visualization plots
   - Important for dosage-sensitive genes (e.g., BRCA1, TP53)

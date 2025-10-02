# GPU-Accelerated Genomic Analysis Pipeline

GPU-accelerated whole-genome analysis pipeline using GPU-native programs for faster processing.

## Pipeline Overview

```
input: FASTQ files
    ↓
1. quality control (FastQC/MultiQC)
    ↓
2. alignment (NVIDIA Parabricks fq2bam) → BAM
    ↓
3. mark duplicates (NVIDIA Parabricks MarkDupliactes)
    ↓
4. base quality score recalibration (NVIDIA Parabricks BQSR)
    ↓
5. variant calling (NVIDIA Parabricks HaplotypeCaller, DeepVariant)
    ↓
6. output VCF!
    ↓
7. variant filtering (NVIDIA Parabricks VariantFiltration)
    ↓
8. variant consequences (SnpEff)
    ↓
9. variant prioritization (bcftools + awk, QUAL>=30)
    ↓
10. structural variant detection (cuteSV + CUDA)
    ↓
11. copy number variant detection (NVIDIA Parabricks GermlineCNVCaller)
    ↓
output: lots of information about your genome!
```

## GPU Requirements

- NVIDIA GPU with CUDA support
- CUDA toolkit installed
- GPU memory: 8GB+ recommended
- GPU compute capability: 6.0+

## Quick Start

### Step 1: Configure for Your Environment
Edit `CONFIG.sh`:
```bash
export SAMPLE_NAME="your_sample_name"
export WORK_DIR="$HOME/wgs_data"
export THREADS=$(nproc)  # or set manually, e.g., 24
export REFERENCE_GENOME="GRCh38"
export GPU_ENABLED=true
```

### Step 2: Run Complete Pipeline
```bash
source CONFIG.sh  # load your configuration
./scripts/gpu/RUN_FULL_PIPELINE_GPU.sh
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
    ├── gpu_sv_caller/
    │   └── {SAMPLE_NAME}.sv.vcf.gz              # Structural variants (DEL/DUP/INV/BND)
    └── cnvkit/
        └── *.cns, *.cnr, *.pdf                  # Copy number variants + plots
```

Note: `{SAMPLE_NAME}` will be replaced with your configured sample name.

## Time Estimates

| Step | Estimated Time |
|------|---------------|
| Quality Control | 30 min |
| Alignment (GPU) | 6-8 hours |
| Mark Duplicates | 1 hour |
| BQSR | 2-3 hours |
| Variant Calling (GPU) | 1-2 hours |
| Filtering | 30 min |
| SnpEff Annotation | 30-45 min |
| Variant Prioritization | 1 min |
| GPU SV Detection | 30-60 min |
| CNVkit CNV | 30 min |
| **TOTAL** | **~12-16 hours** |

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

4. **`gpu_sv_caller/{SAMPLE_NAME}.sv.vcf.gz`**: STRUCTURAL VARIANTS
   - Deletions (DEL), duplications (DUP), inversions (INV), translocations (BND)
   - Check for known disease-associated SVs >50bp

5. **`cnvkit/`**: COPY NUMBER VARIANTS
   - Gene duplications/deletions with visualization plots
   - Important for dosage-sensitive genes (e.g., BRCA1, TP53)

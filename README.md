![Banner](assets/github_banner.png)

# Whole-genome analysis pipeline
This repository contains code to do a full analysis on your whole genome. If you want to run it locally, the repo assumes you have a Linux-based OS, an 8-core CPU with 16 threads, 32GB RAM (64 recommended), 1.5TB of storage (your genome uncompressed will be ~275GB in FASTQ format, even before all the outputs from analysis). This will take about ~30 hours to run the complete end-to-end pipeline. 

If you want to run this on cloud, you just need to follow the steps for running it on cloud. Before you execute, I recommend you get as many free AWS credits as possible. You can get [$100](https://aws.amazon.com/free/offers/) from just signing up with a new credit card + new email, then an additional $100 from completing 5 small challenges on the homescreen (message me if you can't find them). 

---

# Pipeline Overview

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

---

# Analysis Options

## 1. CPU Analysis (`cpu/`)
Complete end-to-end CPU-based genomic analysis pipeline for running on your own Linux workstation.

**Time:** ~21-25 hours  
**Requirements:** Linux workstation with 8+ cores, 32GB+ RAM, 1.5TB+ storage

[View CPU Pipeline Details →](cpu/README.md)

## 2. GPU Analysis (`gpu/`)
GPU-accelerated whole-genome analysis pipeline using GPU-native programs for faster processing.

**Time:** ~12-16 hours  
**Requirements:** NVIDIA GPU with CUDA support, 8GB+ GPU memory

[View GPU Pipeline Details →](gpu/README.md)

## 3. Cloud Infrastructure (`cloud/`)
AWS cloud infrastructure deployment using Terraform for HIPAA/NIST compliant genomic analysis.

**Cost:** ~$300/month  
**Requirements:** AWS account, Terraform installed

[View Cloud Infrastructure Details →](cloud/README.md)

---

## Quick Start

### Step 1: Configure for Your Environment
Edit `CONFIG.sh`:
```bash
export SAMPLE_NAME="your_sample_name"
export WORK_DIR="$HOME/wgs_data"
export THREADS=$(nproc)  # or set manually, e.g., 24
export REFERENCE_GENOME="GRCh38"
export GPU_ENABLED=true  # set to false for CPU-only
```

### Step 2: Choose Your Analysis Method

**For CPU Analysis:**
```bash
source CONFIG.sh
./cpu/scripts/RUN_FULL_PIPELINE_CPU.sh
```

**For GPU Analysis:**
```bash
source CONFIG.sh
./gpu/scripts/RUN_FULL_PIPELINE_GPU.sh
```

**For Cloud Deployment:**
```bash
cd cloud/
terraform init
terraform apply
```

---

## Disclaimer

This software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.

Use at your own risk. The authors assume no responsibility for any damages, data loss, or compliance issues that may arise from the use of this configuration. Always test in a non-production environment first and consult with qualified professionals before deploying to production.

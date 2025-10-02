#!/bin/bash

# S3 upload script
# this script uploads your genomic files to the S3 bucket created by Terraform
# run this from your local machine after downloading your genome files and running Terraform apply

# exit on any error
set -e

# configuration (these will automatically be read from terraform.tfvars)
# you can also set them manually if there are issues
TERRAFORM_DIR="$(dirname "$0")/.."
TERRAFORM_VARS_FILE="$TERRAFORM_DIR/terraform.tfvars"

# read .tfvars file
read_terraform_var() {
    local var_name="$1"
    local default_value="$2"
    
    if [[ -f "$TERRAFORM_VARS_FILE" ]]; then
        # extract values
        local value=$(grep "^${var_name}[[:space:]]*=" "$TERRAFORM_VARS_FILE" | sed 's/.*=[[:space:]]*"\(.*\)".*/\1/' | sed 's/.*=[[:space:]]*\([^"]*\).*/\1/')
        if [[ -n "$value" ]]; then
            echo "$value"
            return
        fi
    fi
    
    # fallback to environment variable or default
    echo "${!var_name:-$default_value}"
}

# read configuration from .tfvars file
BUCKET=$(read_terraform_var "s3_bucket_name" "")
REGION=$(read_terraform_var "region" "")

# validate required variables
if [[ -z "$BUCKET" ]]; then
    echo "Error: S3 bucket name not found in terraform.tfvars or S3_BUCKET_NAME environment variable"
    echo "Please ensure s3_bucket_name is set in your terraform.tfvars file"
    exit 1
fi

# user genome files- customize if this is broken
# ChatGPT can help with this
DEST_PREFIX="${DEST_PREFIX:-wgs_$(date +%Y-%m-%d)}"  # default to today's date
SRC="${SRC:-$HOME/Downloads}"  # default to Downloads folder

# file patterns to upload (customize these for your specific files)
FASTQ_PATTERN="${FASTQ_PATTERN:-*.fastq.gz}"
VCF_PATTERN="${VCF_PATTERN:-*.vcf.gz}"

echo "~configuration:"
echo "  ~bucket: $BUCKET"
echo "  ~region: $REGION"
echo "  ~source: $SRC"
echo "  ~destination: s3://$BUCKET/$DEST_PREFIX/"
echo "  ~FASTQ Pattern: $FASTQ_PATTERN"
echo "  ~VCF Pattern: $VCF_PATTERN"
echo ""

# change to source directory
cd "$SRC"

# check if AWS CLI is installed and configured
if ! command -v aws &> /dev/null; then
    echo "~error: AWS CLI is not installed, please install it first"
    exit 1
fi

# check AWS credentials
if ! aws sts get-caller-identity &> /dev/null; then
    echo "~error: AWS credentials not configured, please run 'aws configure' first."
    exit 1
fi

# try to enable transfer acceleration (optional - will continue if it fails)
echo "~attempting to enable transfer acceleration"
if aws s3api put-bucket-accelerate-configuration \
    --bucket "$BUCKET" \
    --accelerate-configuration Status=Enabled \
    --region "$REGION" 2>/dev/null; then
    echo "~transfer acceleration enabled successfully"
    aws configure set default.s3.use_accelerate_endpoint true
else
    echo "~transfer acceleration not available, continuing with standard endpoint"
    aws configure set default.s3.use_accelerate_endpoint false
fi

# upload files with better timeout and retry settings
echo "~starting upload of genomic files"

# proceed with chunked upload (3 files at a time)
echo "~proceeding with chunked genomic file upload (3 files at a time)"

# create array of files to upload
fastq_files=($(ls "$SRC"/*.fastq.gz 2>/dev/null))
vcf_files=($(ls "$SRC"/*.vcf.gz 2>/dev/null))
all_files=("${fastq_files[@]}" "${vcf_files[@]}")

total_files=${#all_files[@]}
echo "~found $total_files files to upload"

if [[ $total_files -eq 0 ]]; then
    echo "~no genomic files found to upload"
    exit 0
fi

# upload in chunks of 3
chunk_size=3
for ((i=0; i<total_files; i+=chunk_size)); do
    chunk_num=$((i/chunk_size + 1))
    echo ""
    echo "--uploading chunk $chunk_num (files $((i+1))-$((i+chunk_size > total_files ? total_files : i+chunk_size))) of $(( (total_files + chunk_size - 1) / chunk_size ))"
    
    # upload current chunk
    for ((j=i; j<i+chunk_size && j<total_files; j++)); do
        file="${all_files[j]}"
        echo "----uploading: $(basename "$file")"
        aws s3 cp "$file" "s3://$BUCKET/$DEST_PREFIX/" \
            --cli-read-timeout 0 \
            --cli-connect-timeout 300 \
            --region "$REGION"
        
        if [[ $? -eq 0 ]]; then
            echo "  ~✓ uploaded successfully"
        else
            echo "  ~✗ upload failed for $(basename "$file")"
        fi
    done
    
    # pause between chunks (optional)
    if [[ $((i+chunk_size)) -lt $total_files ]]; then
        echo "  ~pausing 2 seconds before next chunk"
        sleep 2
    fi
done

# check what was uploaded
echo ""
echo "~files uploaded to S3:"
aws s3 ls "s3://$BUCKET/$DEST_PREFIX/" --recursive --human-readable --region "$REGION"

echo ""
echo "~upload complete, files are now available in s3://$BUCKET/$DEST_PREFIX/"

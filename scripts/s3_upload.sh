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

# enable transfer acceleration
aws s3api put-bucket-accelerate-configuration \
    --bucket "$BUCKET" \
    --accelerate-configuration Status=Enabled \
    --region "$REGION"

aws s3api get-bucket-accelerate-configuration --bucket "$BUCKET" --region "$REGION"

# configure the AWS CLI for optimal upload performance
aws configure set default.s3.multipart_threshold 64MB
aws configure set default.s3.multipart_chunksize 64MB
aws configure set default.s3.max_concurrent_requests 50
aws configure set default.s3.use_accelerate_endpoint true

# upload files
aws s3 sync "$SRC" "s3://$BUCKET/$DEST_PREFIX/" \
    --exclude "*" \
    --include "$FASTQ_PATTERN" \
    --include "$VCF_PATTERN" \
    --cli-read-timeout 0 \
    --cli-connect-timeout 60 \
    --region "$REGION"

# check what was uploaded
echo ""
echo "~files uploaded to S3:"
aws s3 ls "s3://$BUCKET/$DEST_PREFIX/" --recursive --human-readable --region "$REGION"

echo ""
echo "~upload complete, files are now available in s3://$BUCKET/$DEST_PREFIX/"

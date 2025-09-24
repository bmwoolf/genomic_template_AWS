#!/bin/bash

# your S3 bucket name
S3_BUCKET="s3://your-bucket-name-here/"

# your source directory (where your genetic files are located, usually Downloads)
SOURCE_DIR="${HOME}/Downloads"

# file patterns (customize based on your data, but should be the same)
FASTQ_PATTERN="*.fastq.gz"      # FASTQ files
VCF_PATTERN="*.vcf.gz"          # VCF files

# directory structure in S3
FASTQ_DIR="raw-data/"           # where FASTQ files go
VCF_DIR="results/"              # where VCF files go

echo "uploading genomic data to S3"
echo "bucket: ${S3_BUCKET}"
echo "source: ${SOURCE_DIR}"

# upload FASTQ files
echo "uploading FASTQ files"
aws s3 cp ${SOURCE_DIR}/ ${S3_BUCKET}${FASTQ_DIR} --recursive --include "${FASTQ_PATTERN}" --exclude "*" --cli-read-timeout 0 --cli-write-timeout 0

# upload VCF files
echo "uploading VCF files"
aws s3 cp ${SOURCE_DIR}/ ${S3_BUCKET}${VCF_DIR} --recursive --include "${VCF_PATTERN}" --exclude "*" --cli-read-timeout 0 --cli-write-timeout 0

echo "✅ upload complete!"
echo "Check your S3 bucket: $(echo ${S3_BUCKET} | sed 's|s3://||' | sed 's|/||')"

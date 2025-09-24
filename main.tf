# =============================================================================
# MAIN TERRAFORM CONFIGURATION
# =============================================================================

# local values
locals {
  name = var.project_name
  azs  = ["${var.region}a", "${var.region}b"]
}

# this file now serves as the main entry point
# all infrastructure components have been organized into the components/ folder:
#
# components/
# ├── networking.tf    : VPC, subnets, routing, endpoints
# ├── storage.tf       : S3 buckets, EFS, storage resources
# ├── security.tf      : security groups, KMS keys
# ├── compute.tf       : EC2 instances, user data
# ├── iam.tf           : IAM roles and policies
# ├── monitoring.tf    : CloudWatch, alarms, logging
# ├── compliance.tf    : GuardDuty, Inspector, Config, etc
# └── outputs.tf       : output values
#
# root level files:
# - main.tf            : this file (locals and core resources)
# - variables.tf       : input variables (needed by main.tf)
# - providers.tf       : provider configuration

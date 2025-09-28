# =============================================================================
# MAIN TERRAFORM CONFIGURATION
# =============================================================================

# this file now serves as the main entry point
# local values
locals {
  name = var.project_name
  azs  = ["${var.region}a", "${var.region}b"]
}
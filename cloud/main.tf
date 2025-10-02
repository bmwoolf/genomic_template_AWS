# =============================================================================
# MAIN TERRAFORM CONFIGURATION
# =============================================================================

# local values
locals {
  name = var.project_name
  azs  = ["${var.region}a", "${var.region}b"]
}
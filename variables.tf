variable "region" {
  type    = string
  default = "us-west-2"
}

variable "project_name" {
  type    = string
  default = "wgs-paid"
}

variable "home_ip_cidr" {
  type        = string
  description = "Your home IP CIDR block (e.g. 203.0.113.42/32)"
}

variable "ssh_key_name" {
  type        = string
  description = "Existing EC2 key pair name"
}

variable "s3_bucket_name" {
  type        = string
  description = "Globally unique S3 bucket name"
}

variable "instance_type" {
  type    = string
  default = "c7i.4xlarge" # pick one: c7i.4xlarge | c7i.8xlarge | g5.2xlarge
  validation {
    condition     = contains(["c7i.4xlarge", "g5.2xlarge"], var.instance_type)
    error_message = "instance_type must be one of: c7i.4xlarge, g5.2xlarge."
  }
}

variable "disk_gb" {
  type        = number
  default     = 500
  description = "Scratch space for BAM/VCF files"
}

variable "use_spot" {
  type        = bool
  default     = true
  description = "Use spot instances for cost savings"
}

# OPTIONAL: if true, a sudo shutdown from the OS will TERMINATE the instance,
# which also deletes the root/scratch EBS. If false, shutdown just stops it.
variable "terminate_on_shutdown" {
  type    = bool
  default = false
}

variable "alert_email" {
  type        = string
  description = "Email address for CloudWatch and budget alerts"
}

# Sensitive variables for secrets and passwords
variable "database_password" {
  type        = string
  description = "Database password for the application"
  sensitive   = true
}


variable "application_secret" {
  type        = string
  description = "Application secret key for authentication"
  sensitive   = true
}

variable "jwt_secret" {
  type        = string
  description = "JWT secret key for token signing"
  sensitive   = true
  default     = ""
}

variable "encryption_key" {
  type        = string
  description = "Application-level encryption key"
  sensitive   = true
  default     = ""
}
# =============================================================================
# COMPLIANCE AND SECURITY SERVICES
# =============================================================================

# CloudTrail for audit logging
resource "aws_cloudtrail" "main" {
  name                          = "${local.name}-cloudtrail"
  s3_bucket_name                = aws_s3_bucket.cloudtrail.id
  kms_key_id                    = aws_kms_key.cloudtrail.arn
  include_global_service_events = true
  is_multi_region_trail         = true
  enable_logging                = true
  enable_log_file_validation    = true

  event_selector {
    read_write_type           = "All"
    include_management_events = true

    data_resource {
      type   = "AWS::S3::Object"
      values = ["arn:aws:s3:::${aws_s3_bucket.fastq.bucket}/"]
    }
  }

  tags = { Name = "${local.name}-cloudtrail" }
}

# AWS Config for compliance
resource "aws_config_configuration_recorder" "main" {
  name     = "${local.name}-config-recorder"
  role_arn = aws_iam_role.config_role.arn

  recording_group {
    all_supported                 = true
    include_global_resource_types = true
  }
}

resource "aws_config_delivery_channel" "main" {
  name           = "${local.name}-config-delivery"
  s3_bucket_name = aws_s3_bucket.config.id
  s3_key_prefix  = "config"
}

resource "aws_config_configuration_recorder_status" "main" {
  name       = aws_config_configuration_recorder.main.name
  is_enabled = true
  depends_on = [aws_config_delivery_channel.main]
}

# AWS Config Rules for compliance monitoring
resource "aws_config_config_rule" "s3_encryption" {
  name = "s3-bucket-server-side-encryption-enabled"
  source {
    owner             = "AWS"
    source_identifier = "S3_BUCKET_SERVER_SIDE_ENCRYPTION_ENABLED"
  }
  depends_on = [aws_config_configuration_recorder.main]
}

resource "aws_config_config_rule" "imdsv2" {
  name = "ec2-imdsv2-enabled"
  source {
    owner             = "AWS"
    source_identifier = "EC2_IMDSV2_CHECK"
  }
  depends_on = [aws_config_configuration_recorder.main]
}

resource "aws_config_config_rule" "s3_public_access" {
  name = "s3-bucket-public-access-prohibited"
  source {
    owner             = "AWS"
    source_identifier = "S3_BUCKET_PUBLIC_ACCESS_PROHIBITED"
  }
  depends_on = [aws_config_configuration_recorder.main]
}

resource "aws_config_config_rule" "ebs_encryption" {
  name = "ebs-encryption-enabled"
  source {
    owner             = "AWS"
    source_identifier = "EBS_ENCRYPTION_ENABLED"
  }
  depends_on = [aws_config_configuration_recorder.main]
}

resource "aws_config_config_rule" "cloudtrail_enabled" {
  name = "cloudtrail-enabled"
  source {
    owner             = "AWS"
    source_identifier = "CLOUD_TRAIL_ENABLED"
  }
  depends_on = [aws_config_configuration_recorder.main]
}

# AWS GuardDuty for threat detection
resource "aws_guardduty_detector" "main" {
  enable                       = true
  finding_publishing_frequency = "FIFTEEN_MINUTES"
  tags                         = { Name = "${local.name}-guardduty" }
}

# AWS Inspector for vulnerability assessment
resource "aws_inspector2_enabler" "main" {
  account_ids    = [data.aws_caller_identity.current.account_id]
  resource_types = ["EC2", "ECR", "LAMBDA"]
}

# AWS Inspector assessment target
resource "aws_inspector_assessment_target" "main" {
  name               = "${local.name}-assessment-target"
  resource_group_arn = aws_inspector_resource_group.main.arn
}

# AWS Inspector resource group
resource "aws_inspector_resource_group" "main" {
  tags = {
    Name = "${local.name}-inspector-resources"
  }
}

# AWS Inspector assessment template
resource "aws_inspector_assessment_template" "main" {
  name       = "${local.name}-assessment-template"
  target_arn = aws_inspector_assessment_target.main.arn
  duration   = 3600

  rules_package_arns = [
    "arn:aws:inspector:${var.region}:316112463485:rulespackage/0-9hgA516P",
    "arn:aws:inspector:${var.region}:316112463485:rulespackage/0-H5hpSawc",
    "arn:aws:inspector:${var.region}:316112463485:rulespackage/0-JJOtZiqQ",
    "arn:aws:inspector:${var.region}:316112463485:rulespackage/0-cE5K8AlX"
  ]
}

# AWS Macie for data discovery and protection
resource "aws_macie2_account" "main" {
  finding_publishing_frequency = "FIFTEEN_MINUTES"
  status                       = "ENABLED"
}

# AWS Macie classification job for S3 bucket
resource "aws_macie2_classification_job" "fastq_bucket" {
  job_type = "ONE_TIME"
  name     = "${local.name}-fastq-classification-job"

  s3_job_definition {
    bucket_definitions {
      account_id = data.aws_caller_identity.current.account_id
      buckets    = [aws_s3_bucket.fastq.bucket]
    }
  }

  depends_on = [aws_macie2_account.main]
}

# AWS Security Hub for centralized security findings
resource "aws_securityhub_account" "main" {
  enable_default_standards = true
}

# AWS Security Hub standards subscriptions
resource "aws_securityhub_standards_subscription" "cis" {
  standards_arn = "arn:aws:securityhub:${var.region}::standards/cis-aws-foundations-benchmark/v/1.2.0"
  depends_on    = [aws_securityhub_account.main]
}

resource "aws_securityhub_standards_subscription" "pci" {
  standards_arn = "arn:aws:securityhub:${var.region}::standards/pci-dss/v/3.2.1"
  depends_on    = [aws_securityhub_account.main]
}

# AWS Systems Manager Session Manager for secure shell access
resource "aws_ssm_document" "session_manager_preferences" {
  name            = "SSM-SessionManagerRunShell"
  document_type   = "Session"
  document_format = "JSON"

  content = jsonencode({
    schemaVersion = "1.0"
    description   = "Document to hold regional settings for Session Manager"
    sessionType   = "Standard_Stream"
    inputs = {
      s3BucketName                = aws_s3_bucket.cloudtrail.bucket
      s3KeyPrefix                 = "session-manager"
      s3EncryptionEnabled         = true
      cloudWatchLogGroupName      = aws_cloudwatch_log_group.ec2_logs.name
      cloudWatchEncryptionEnabled = true
      kmsKeyId                    = aws_kms_key.cloudtrail.arn
    }
  })
}

# AWS Systems Manager Parameter Store for secure configuration
resource "aws_ssm_parameter" "database_password" {
  name   = "/${local.name}/database/password"
  type   = "SecureString"
  value  = var.database_password
  key_id = aws_kms_key.cloudtrail.key_id

  tags = {
    Name       = "${local.name}-db-password"
    Compliance = "HIPAA"
  }
}

# AWS Backup for data protection and recovery
resource "aws_backup_vault" "main" {
  name        = "${local.name}-backup-vault"
  kms_key_arn = aws_kms_key.cloudtrail.arn
  tags        = { Name = "${local.name}-backup-vault" }
}

# AWS Backup plan
resource "aws_backup_plan" "main" {
  name = "${local.name}-backup-plan"

  rule {
    rule_name         = "daily_backup"
    target_vault_name = aws_backup_vault.main.name
    schedule          = "cron(0 2 * * ? *)" # Daily at 2 AM

    lifecycle {
      cold_storage_after = 30
      delete_after       = 120
    }

    recovery_point_tags = {
      Name = "${local.name}-backup"
    }
  }

  tags = { Name = "${local.name}-backup-plan" }
}

# AWS Backup selection for EC2 instances
resource "aws_backup_selection" "ec2" {
  iam_role_arn = aws_iam_role.backup_role.arn
  name         = "${local.name}-ec2-backup-selection"
  plan_id      = aws_backup_plan.main.id

  resources = [
    aws_instance.worker.arn,
    aws_instance.bastion.arn
  ]
}

# AWS WAF for web application protection
resource "aws_wafv2_web_acl" "main" {
  name  = "${local.name}-waf"
  scope = "REGIONAL"

  default_action {
    allow {}
  }

  rule {
    name     = "AWSManagedRulesCommonRuleSet"
    priority = 1

    override_action {
      none {}
    }

    statement {
      managed_rule_group_statement {
        name        = "AWSManagedRulesCommonRuleSet"
        vendor_name = "AWS"
      }
    }

    visibility_config {
      cloudwatch_metrics_enabled = true
      metric_name                = "CommonRuleSetMetric"
      sampled_requests_enabled   = true
    }
  }

  visibility_config {
    cloudwatch_metrics_enabled = true
    metric_name                = "${local.name}WAFMetric"
    sampled_requests_enabled   = true
  }

  tags = { Name = "${local.name}-waf" }
}

# AWS Certificate Manager for SSL/TLS certificates
resource "aws_acm_certificate" "main" {
  domain_name       = "*.${local.name}.com"
  validation_method = "DNS"

  lifecycle {
    create_before_destroy = true
  }

  tags = { Name = "${local.name}-certificate" }
}

# AWS Secrets Manager for sensitive data
resource "aws_secretsmanager_secret" "main" {
  name                    = "${local.name}-secrets"
  description             = "Secrets for ${local.name} application"
  kms_key_id              = aws_kms_key.cloudtrail.key_id
  recovery_window_in_days = 30

  tags = {
    Name       = "${local.name}-secrets"
    Compliance = "HIPAA"
  }
}

resource "aws_secretsmanager_secret_version" "main" {
  secret_id = aws_secretsmanager_secret.main.id
  secret_string = jsonencode({
    database_password  = var.database_password
    api_key            = var.api_key
    application_secret = var.application_secret
    jwt_secret         = var.jwt_secret != "" ? var.jwt_secret : null
    encryption_key     = var.encryption_key != "" ? var.encryption_key : null
  })
}

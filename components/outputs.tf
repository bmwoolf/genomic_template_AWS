output "bastion_public_ip" { 
  description = "Public IP of the bastion host"
  value       = aws_instance.bastion.public_ip 
}

output "worker_private_ip" { 
  description = "Private IP of the worker instance"
  value       = aws_instance.worker.private_ip 
}

output "ssh_command_bastion" { 
  description = "SSH command to connect to bastion host"
  value       = "ssh -i ~/.ssh/<your-key>.pem ubuntu@${aws_instance.bastion.public_ip}" 
}

output "ssh_command_worker" { 
  description = "SSH command to connect to worker instance via bastion"
  value       = "ssh -i ~/.ssh/<your-key>.pem ubuntu@${aws_instance.worker.private_ip} (via bastion)" 
}

output "bucket" { 
  description = "S3 bucket name for FASTQ files"
  value       = aws_s3_bucket.fastq.bucket 
}

output "efs_id" { 
  description = "EFS file system ID for reference genomes"
  value       = aws_efs_file_system.efs.id 
}

output "cloudtrail_bucket" { 
  description = "S3 bucket name for CloudTrail logs"
  value       = aws_s3_bucket.cloudtrail.bucket 
}

output "config_bucket" { 
  description = "S3 bucket name for AWS Config logs"
  value       = aws_s3_bucket.config.bucket 
}

# Compliance and Security Outputs
output "guardduty_detector_id" {
  description = "GuardDuty detector ID for threat detection"
  value       = aws_guardduty_detector.main.id
}

output "security_hub_arn" {
  description = "Security Hub ARN for centralized security findings"
  value       = aws_securityhub_account.main.arn
}

output "backup_vault_arn" {
  description = "AWS Backup vault ARN for data protection"
  value       = aws_backup_vault.main.arn
}

output "secrets_manager_arn" {
  description = "Secrets Manager ARN for sensitive data storage"
  value       = aws_secretsmanager_secret.main.arn
}

output "kms_keys" {
  description = "KMS key ARNs for encryption"
  value = {
    s3_fastq    = aws_kms_key.s3_fastq.arn
    cloudtrail  = aws_kms_key.cloudtrail.arn
    config      = aws_kms_key.config.arn
    sns         = aws_kms_key.sns.arn
  }
}

output "compliance_status" {
  description = "Compliance services status"
  value = {
    guardduty_enabled     = aws_guardduty_detector.main.enable
    inspector_enabled     = aws_inspector2_enabler.main.resource_types
    macie_enabled         = aws_macie2_account.main.status
    security_hub_enabled  = aws_securityhub_account.main.enable_default_standards
    backup_enabled        = aws_backup_plan.main.name
    config_rules_enabled  = length(aws_config_config_rule.s3_encryption)
    vpc_endpoints_enabled = length(aws_vpc_endpoint.ifaces)
  }
}

output "access_logs_bucket" {
  description = "S3 bucket name for immutable access logs"
  value       = aws_s3_bucket.logs.bucket
}

output "session_manager_info" {
  description = "Session Manager connection information"
  value = "Use AWS CLI: aws ssm start-session --target ${aws_instance.worker.id}"
}

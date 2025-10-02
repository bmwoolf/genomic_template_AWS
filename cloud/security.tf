# =============================================================================
# SECURITY INFRASTRUCTURE
# =============================================================================

# Security groups
resource "aws_security_group" "bastion_sg" {
  name        = "${local.name}-bastion-sg"
  description = "Bastion host SSH access"
  vpc_id      = aws_vpc.vpc.id
  ingress {
    description = "SSH from home"
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = [var.home_ip_cidr]
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
  tags = { Name = "${local.name}-bastion-sg" }
}

resource "aws_security_group" "ec2_sg" {
  name        = "${local.name}-ec2-sg"
  description = "EC2 compute instances"
  vpc_id      = aws_vpc.vpc.id
  ingress {
    description     = "SSH from bastion"
    from_port       = 22
    to_port         = 22
    protocol        = "tcp"
    security_groups = [aws_security_group.bastion_sg.id]
  }
  ingress {
    description = "Custom application ports"
    from_port   = 8000
    to_port     = 8999
    protocol    = "tcp"
    cidr_blocks = ["10.42.0.0/16"]  # Allow from VPC CIDR
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
  tags = { Name = "${local.name}-ec2-sg" }
}

resource "aws_security_group" "efs_sg" {
  name        = "${local.name}-efs-sg"
  description = "EFS"
  vpc_id      = aws_vpc.vpc.id
  ingress {
    from_port       = 2049
    to_port         = 2049
    protocol        = "tcp"
    security_groups = [aws_security_group.ec2_sg.id]
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

# KMS Keys for encryption
resource "aws_kms_key" "s3_fastq" {
  description             = "CMK for FASTQ bucket"
  enable_key_rotation     = true
  deletion_window_in_days = 30
}

resource "aws_kms_alias" "s3_fastq" {
  name          = "alias/${local.name}-s3-fastq"
  target_key_id = aws_kms_key.s3_fastq.key_id
}

# KMS key for CloudTrail
resource "aws_kms_key" "cloudtrail" {
  description             = "CMK for CloudTrail"
  enable_key_rotation     = true
  deletion_window_in_days = 30
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "Enable IAM User Permissions"
        Effect = "Allow"
        Principal = {
          AWS = "arn:aws:iam::${data.aws_caller_identity.current.account_id}:root"
        }
        Action   = "kms:*"
        Resource = "*"
      },
      {
        Sid    = "Allow CloudTrail to encrypt logs"
        Effect = "Allow"
        Principal = {
          Service = "cloudtrail.amazonaws.com"
        }
        Action = [
          "kms:GenerateDataKey*"
        ]
        Resource = "*"
        Condition = {
          StringEquals = {
            "kms:EncryptionContext:aws:cloudtrail:arn" = "arn:aws:cloudtrail:${var.region}:${data.aws_caller_identity.current.account_id}:trail/${local.name}-cloudtrail"
          }
        }
      },
      {
        Sid    = "Allow key policy updates"
        Effect = "Allow"
        Principal = {
          AWS = "arn:aws:iam::${data.aws_caller_identity.current.account_id}:root"
        }
        Action = [
          "kms:PutKeyPolicy"
        ]
        Resource = "*"
      }
    ]
  })
  tags = { Name = "${local.name}-cloudtrail-kms" }
}

resource "aws_kms_alias" "cloudtrail" {
  name          = "alias/${local.name}-cloudtrail"
  target_key_id = aws_kms_key.cloudtrail.key_id
}

# KMS key for AWS Config
resource "aws_kms_key" "config" {
  description             = "CMK for AWS Config"
  enable_key_rotation     = true
  deletion_window_in_days = 30
  tags = { Name = "${local.name}-config-kms" }
}

resource "aws_kms_alias" "config" {
  name          = "alias/${local.name}-config"
  target_key_id = aws_kms_key.config.key_id
}

# KMS key for SNS encryption
resource "aws_kms_key" "sns" {
  description             = "CMK for SNS topic encryption"
  enable_key_rotation     = true
  deletion_window_in_days = 30
  tags = { Name = "${local.name}-sns-kms" }
}

resource "aws_kms_alias" "sns" {
  name          = "alias/${local.name}-sns"
  target_key_id = aws_kms_key.sns.key_id
}

# data source for current AWS account
data "aws_caller_identity" "current" {}

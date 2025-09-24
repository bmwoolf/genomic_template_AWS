# =============================================================================
# COMPUTE INFRASTRUCTURE
# =============================================================================

# AMI (Ubuntu 22.04 LTS)
data "aws_ami" "ubuntu" {
  most_recent = true
  owners      = ["099720109477"] # Canonical
  filter {
    name   = "name"
    values = ["ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server-*"]
  }
}

# User data (Docker + AWS CLI + EFS; GPU setup if present)
locals {
  user_data = <<-CLOUD
    #!/bin/bash
    set -eux
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y docker.io awscli amazon-efs-utils tmux build-essential pciutils unzip

    usermod -aG docker ubuntu
    mkdir -p /data/fastq /data/results /mnt/efs /opt/aws
    chown -R ubuntu:ubuntu /data

    # Install CloudWatch agent
    wget https://s3.amazonaws.com/amazoncloudwatch-agent/ubuntu/amd64/latest/amazon-cloudwatch-agent.deb
    dpkg -i -E ./amazon-cloudwatch-agent.deb

    # Install Systems Manager agent for Session Manager
    wget https://s3.amazonaws.com/ec2-downloads-windows/SSMAgent/latest/debian_amd64/amazon-ssm-agent.deb
    dpkg -i amazon-ssm-agent.deb
    systemctl enable amazon-ssm-agent
    systemctl start amazon-ssm-agent

    # Mount EFS
    echo "${aws_efs_file_system.efs.id}:/ /mnt/efs efs _netdev,tls 0 0" >> /etc/fstab || true
    systemctl daemon-reload || true
    mount -a || true
    chown -R ubuntu:ubuntu /mnt/efs || true

    # If GPU (e.g., g5.2xlarge), install NVIDIA driver + nvidia-docker
    if lspci | grep -i -q 'nvidia'; then
      apt-get install -y ubuntu-drivers-common curl
      ubuntu-drivers autoinstall || true
      # nvidia-docker2
      distribution=$(. /etc/os-release; echo $ID$VERSION_ID) && \
      curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg && \
      curl -fsSL https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list | \
        sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
        tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
      apt-get update
      apt-get install -y nvidia-container-toolkit
      nvidia-ctk runtime configure --runtime=docker
      systemctl restart docker
    fi

    # Pre-pull common images
    docker pull biocontainers/bwa:v0.7.17_cv1 || true
    docker pull biocontainers/samtools:v1.17-1-deb_cv1 || true
    docker pull broadinstitute/gatk:4.5.0.0 || true

    # Start CloudWatch agent
    /opt/aws/amazon-cloudwatch-agent/bin/amazon-cloudwatch-agent-ctl \
      -a fetch-config \
      -m ec2 \
      -c default \
      -s

    # Configure automatic security updates for compliance
    echo 'Unattended-Upgrade::Automatic-Reboot "true";' >> /etc/apt/apt.conf.d/50unattended-upgrades
    echo 'Unattended-Upgrade::Automatic-Reboot-Time "02:00";' >> /etc/apt/apt.conf.d/50unattended-upgrades
    systemctl enable unattended-upgrades
    systemctl start unattended-upgrades
  CLOUD
}

# Bastion host for SSH access to private instances
resource "aws_instance" "bastion" {
  ami                    = data.aws_ami.ubuntu.id
  instance_type          = "t3.micro"
  subnet_id              = aws_subnet.public_a.id
  vpc_security_group_ids = [aws_security_group.bastion_sg.id]
  key_name               = var.ssh_key_name

  root_block_device {
    volume_size           = 20
    volume_type           = "gp3"
    encrypted             = true
    delete_on_termination = true
  }

  tags = { Name = "${local.name}-bastion" }
}

# EC2 (Spot optional)
resource "aws_instance" "worker" {
  ami                    = data.aws_ami.ubuntu.id
  instance_type          = var.instance_type
  subnet_id              = aws_subnet.private_a.id
  vpc_security_group_ids = [aws_security_group.ec2_sg.id]
  key_name               = var.ssh_key_name
  iam_instance_profile   = aws_iam_instance_profile.ec2_profile.name

  # if you want OS shutdown to fully terminate (and auto-delete EBS), flip the var to true
  instance_initiated_shutdown_behavior = var.terminate_on_shutdown ? "terminate" : "stop"

  # be explicit that the root/scratch EBS is deleted when the instance is terminated
  root_block_device {
    volume_size           = var.disk_gb
    volume_type           = "gp3"
    encrypted             = true
    delete_on_termination = true
  }

  # EC2 hardening: IMDSv2 enforcement
  metadata_options {
    http_tokens                 = "required" # enforce IMDSv2
    http_put_response_hop_limit = 1          # prevent container access to metadata
  }

  user_data = local.user_data

  dynamic "instance_market_options" {
    for_each = var.use_spot ? [1] : []
    content {
      market_type = "spot"
      spot_options {
        spot_instance_type             = "one-time"
        instance_interruption_behavior = "terminate"
      }
    }
  }

  tags = {
    Name = "${local.name}-ec2"
    Type = "compute"
  }
}

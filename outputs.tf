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

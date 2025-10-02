# AWS Cloud Infrastructure in Terraform (HIPAA/NIST Compliant)

This Terraform configuration deploys a comprehensive, compliance-ready AWS infrastructure for genetic analysis workloads, meeting HIPAA, NIH-GDS, GINA, and NIST SP 800-171 requirements.

All development is done on a Mac. You may need to change commands depending on your operating system.


## Quick start (Cloud)

1. install the AWS CLI
   ```bash
   brew install aws
   aws configure
   # enter in ACCESS_KEY_ID and SECRET_ACCESS_KEY from AWS console
   # if you don't know how to do that, ChatGPT can help
   ```

2. create an EC2 key pair in the AWS console (ssh_key_name)
![key-pair](../assets/ec2_key_pair.png)

3. configure variables for your environment:
   ```bash
   # copy + edit terraform.tfvars with your values
   cp terraform.tfvars.example terraform.tfvars

   # create database_password (save for next step)
   openssl rand -base64 32

   # create application_secret (save for next step)
   openssl rand -base64 32
   ```

4. fill out the tfvars file with your info
   ```
   region              = "us-east-1"                # use whatever is closest to you
   project_name        = "your-project-name"        # like "wgs-analysis-your-name"
   home_ip_cidr        = "YOUR.IP.ADDRESS/32"       # your public IP with /32
   ssh_key_name        = "your-ec2-keypair-name"
   s3_bucket_name      = "your-unique-bucket-name"  # must be globally unique
   instance_type       = "c7i.4xlarge"              # or "g5.2xlarge" for GPU
   disk_gb             = 500
   use_spot            = true
   alert_email         = "your-email@example.com"
   database_password   = "your-secure-database-password-here"
   application_secret  = "your-application-secret-key-for-auth"
   ```

5. deploy infrastructure:
   ```bash
   terraform fmt                    # makes it pretty
   terraform init                   # connects to your provider (AWS)
   terraform plan -out=plan.tfplan  # creates and saves the plan
   terraform apply                  # deploys infra

   # take down infra
   terraform destroy                # takes down infra
   ```

6. access EC2 instances to run software on them (for your genome):
   ```bash
   # via Session Manager (recommended, more compliant)
   aws ssm start-session --target <instance-id>
   
   # via Bastion Host (traditional)
   ssh -i ~/.ssh/<key>.pem ubuntu@<bastion-ip>
   ```

7. upload genome to S3 via S3 Transfer Acceleration and parallelism
   ```bash
   ./scripts/s3_upload.sh
   ```


### Cost estimation
- monthly cost: ~$300
- major components: EC2 instances, NAT Gateway, CloudWatch logs, compliance services
- cost optimization: Spot instances, S3 lifecycle policies, log retention

### Resources
- total resources: 110
- VPC endpoints: 10 (complete private connectivity)
- S3 buckets: 4 (FASTQ, logs, CloudTrail, Config)
- KMS keys: 4 (dedicated encryption per service)
- Config rules: 5 (continuous compliance monitoring)
- Security services: 6 (GuardDuty, Inspector, Macie, Security Hub, WAF, Backup)


###  Files
- `main.tf`: main entry point with locals and documentation
- `variables.tf`: input variables and validation
- `providers.tf`: AWS provider configuration
- `terraform.tfvars.example`: configuration template
- `README.md`: this documentation
- `networking.tf`: VPC, subnets, routing, NAT Gateway, VPC endpoints
- `storage.tf`: S3 buckets (FASTQ, logs, CloudTrail, Config), EFS file system
- `security.tf`: security groups, KMS keys, encryption policies
- `compute.tf`: EC2 instances, AMI selection, user data scripts
- `iam.tf`: IAM roles, policies, instance profiles
- `monitoring.tf`: CloudWatch logs, alarms, SNS topics, budgets
- `compliance.tf`: GuardDuty, Inspector, Config, Macie, Security Hub, Backup
- `outputs.tf`: output values and resource information

---

## Compliance for HIPAA, NIH-GDS, GINA, and NIST SP 800-171

To comply with HIPAA, NIH-GDS, GINA, and NIST SP 800-171 for genetic analysis workloads, you must encrypt all data at rest (S3, EBS, EFS using KMS) and in transit (TLS enforcement via VPC endpoints), maintain immutable audit logs with object lock for 365 days, implement least-privilege IAM policies, isolate compute in private subnets with NAT Gateway, enable continuous security monitoring through GuardDuty, Inspector, Macie, Security Hub, and WAF, configure multi-region CloudTrail logging with automated Config compliance checks, capture VPC flow logs, harden instance metadata to IMDSv2, use Session Manager for SSH-less access, and establish automated daily backups- all of which this Terraform configuration implements by default, though you should review with your organization's compliance team and legal counsel to ensure it meets your specific regulatory requirements. Note that it is not cheap- I was running fully-compliant infrastructure for this single genome and it costs ~$10/day.


---

## Compliance links
- AWS compliance programs: https://aws.amazon.com/compliance/
- HIPAA on AWS: https://aws.amazon.com/compliance/hipaa-compliance/
- NIST SP 800-171: https://csrc.nist.gov/publications/detail/sp/800-171/rev-2/final
- NIH-GDS policy: https://gds.nih.gov/

---

Note: This configuration provides enterprise-grade security and compliance. Please review with your organization's compliance team and legal counsel to ensure it meets all specific requirements for your use case.

## Disclaimer

This software is provided "as is" without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.

Use at your own risk. The authors assume no responsibility for any damages, data loss, or compliance issues that may arise from the use of this configuration. Always test in a non-production environment first and consult with qualified professionals before deploying to production.

# genomic_template_AWS
A Terraform repo that creates the necessary infrastructure in AWS to analyze your own genome.

1. install the AWS CLI (Homebrew: `brew install aws`)
2. configure your AWS CLI locally via `aws configure`
3. fill out the `terraform.tfvars` with your info to deploy the repo

The goal is to be HIPAA, NIH-GDS, GINA, and NIST SP 800-171 compliant. It wil probably need geographical tweaks if you are from the EU so you can abide by GDPR.
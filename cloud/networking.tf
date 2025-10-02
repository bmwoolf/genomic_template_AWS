# =============================================================================
# NETWORKING INFRASTRUCTURE
# =============================================================================

# VPC
resource "aws_vpc" "vpc" {
  cidr_block           = "10.42.0.0/16"
  enable_dns_hostnames = true
  tags = { Name = "${local.name}-vpc" }
}

# internet gateway
resource "aws_internet_gateway" "igw" {
  vpc_id = aws_vpc.vpc.id
  tags = { Name = "${local.name}-igw" }
}

# public subnets for NAT gateway and bastion host
resource "aws_subnet" "public_a" {
  vpc_id                  = aws_vpc.vpc.id
  cidr_block              = "10.42.1.0/24"
  availability_zone       = local.azs[0]
  map_public_ip_on_launch = true
  tags = { 
    Name = "${local.name}-public-a"
    Type = "public"
  }
}

resource "aws_subnet" "public_b" {
  vpc_id                  = aws_vpc.vpc.id
  cidr_block              = "10.42.2.0/24"
  availability_zone       = local.azs[1]
  map_public_ip_on_launch = true
  tags = { 
    Name = "${local.name}-public-b"
    Type = "public"
  }
}

# private subnets for compute resources
resource "aws_subnet" "private_a" {
  vpc_id            = aws_vpc.vpc.id
  cidr_block        = "10.42.11.0/24"
  availability_zone = local.azs[0]
  tags = { 
    Name = "${local.name}-private-a"
    Type = "private"
  }
}

resource "aws_subnet" "private_b" {
  vpc_id            = aws_vpc.vpc.id
  cidr_block        = "10.42.12.0/24"
  availability_zone = local.azs[1]
  tags = { 
    Name = "${local.name}-private-b"
    Type = "private"
  }
}

# route tables
resource "aws_route_table" "public" {
  vpc_id = aws_vpc.vpc.id
  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.igw.id
  }
  tags = { Name = "${local.name}-public-rt" }
}

resource "aws_route_table" "private" {
  vpc_id = aws_vpc.vpc.id
  # No internet route - private subnet uses bastion proxy for outbound access
  tags = { Name = "${local.name}-private-rt" }
}

# route table associations
resource "aws_route_table_association" "public_a" {
  subnet_id      = aws_subnet.public_a.id
  route_table_id = aws_route_table.public.id
}

resource "aws_route_table_association" "public_b" {
  subnet_id      = aws_subnet.public_b.id
  route_table_id = aws_route_table.public.id
}

resource "aws_route_table_association" "private_a" {
  subnet_id      = aws_subnet.private_a.id
  route_table_id = aws_route_table.private.id
}

resource "aws_route_table_association" "private_b" {
  subnet_id      = aws_subnet.private_b.id
  route_table_id = aws_route_table.private.id
}

# NAT Gateway removed for cost savings - using bastion proxy instead

# VPC endpoints
# S3 gateway VPC endpoint (keep traffic on AWS backbone)
resource "aws_vpc_endpoint" "s3_gw" {
  vpc_id            = aws_vpc.vpc.id
  service_name      = "com.amazonaws.${var.region}.s3"
  vpc_endpoint_type = "Gateway"
  route_table_ids   = [aws_route_table.public.id, aws_route_table.private.id]
  tags = { Name = "${local.name}-s3-endpoint" }
}

# VPC interface endpoints for private AWS service access
# only keep essential endpoints for genomic analysis
locals {
  endpoints = [
    "com.amazonaws.${var.region}.ecr.api",    # Container registry API
    "com.amazonaws.${var.region}.ecr.dkr",   # Container registry Docker
    "com.amazonaws.${var.region}.kms"        # Encryption keys
  ]
}

resource "aws_vpc_endpoint" "ifaces" {
  for_each = toset(local.endpoints)
  
  vpc_id               = aws_vpc.vpc.id
  service_name         = each.value
  vpc_endpoint_type    = "Interface"
  subnet_ids           = [aws_subnet.private_a.id, aws_subnet.private_b.id]
  security_group_ids   = [aws_security_group.ec2_sg.id]
  private_dns_enabled  = true
  
  tags = { Name = "${local.name}-${replace(each.value, "com.amazonaws.${var.region}.", "")}" }
}

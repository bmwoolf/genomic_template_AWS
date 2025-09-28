# =============================================================================
# MONITORING AND LOGGING INFRASTRUCTURE
# =============================================================================

# CloudWatch log groups
resource "aws_cloudwatch_log_group" "ec2_logs" {
  name              = "/aws/ec2/${local.name}"
  retention_in_days = 30
  tags              = { Name = "${local.name}-ec2-logs" }
}

resource "aws_cloudwatch_log_group" "vpc_flow_logs" {
  name              = "/aws/vpc/flowlogs/${local.name}"
  retention_in_days = 14
  tags              = { Name = "${local.name}-vpc-flow-logs" }
}

# VPC flow logs
resource "aws_flow_log" "vpc_flow_logs" {
  iam_role_arn    = aws_iam_role.vpc_flow_logs_role.arn
  log_destination = aws_cloudwatch_log_group.vpc_flow_logs.arn
  traffic_type    = "ALL"
  vpc_id          = aws_vpc.vpc.id
  tags            = { Name = "${local.name}-vpc-flow-logs" }
}

# CloudWatch alarms
resource "aws_cloudwatch_metric_alarm" "high_cpu" {
  alarm_name          = "${local.name}-high-cpu"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "CPUUtilization"
  namespace           = "AWS/EC2"
  period              = "300"
  statistic           = "Average"
  threshold           = "80"
  alarm_description   = "This metric monitors ec2 cpu utilization"
  alarm_actions       = [aws_sns_topic.alerts.arn]
  dimensions = {
    InstanceId = aws_instance.worker.id
  }
  tags = { Name = "${local.name}-high-cpu-alarm" }
}

resource "aws_cloudwatch_metric_alarm" "high_memory" {
  alarm_name          = "${local.name}-high-memory"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "MemoryUtilization"
  namespace           = "CWAgent"
  period              = "300"
  statistic           = "Average"
  threshold           = "85"
  alarm_description   = "This metric monitors ec2 memory utilization"
  alarm_actions       = [aws_sns_topic.alerts.arn]
  dimensions = {
    InstanceId = aws_instance.worker.id
  }
  tags = { Name = "${local.name}-high-memory-alarm" }
}

resource "aws_cloudwatch_metric_alarm" "high_disk" {
  alarm_name          = "${local.name}-high-disk"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "disk_used_percent"
  namespace           = "CWAgent"
  period              = "300"
  statistic           = "Average"
  threshold           = "90"
  alarm_description   = "This metric monitors ec2 disk utilization"
  alarm_actions       = [aws_sns_topic.alerts.arn]
  dimensions = {
    InstanceId = aws_instance.worker.id
    device     = "/dev/nvme0n1"
    fstype     = "ext4"
    path       = "/"
  }
  tags = { Name = "${local.name}-high-disk-alarm" }
}

# GPU monitoring (if GPU instance)
resource "aws_cloudwatch_metric_alarm" "gpu_utilization" {
  count               = var.instance_type == "g5.2xlarge" ? 1 : 0
  alarm_name          = "${local.name}-gpu-utilization"
  comparison_operator = "GreaterThanThreshold"
  evaluation_periods  = "2"
  metric_name         = "GPUUtilization"
  namespace           = "CWAgent"
  period              = "300"
  statistic           = "Average"
  threshold           = "95"
  alarm_description   = "This metric monitors GPU utilization"
  alarm_actions       = [aws_sns_topic.alerts.arn]
  dimensions = {
    InstanceId = aws_instance.worker.id
  }
  tags = { Name = "${local.name}-gpu-utilization-alarm" }
}

# SNS topic for alerts with KMS encryption
resource "aws_sns_topic" "alerts" {
  name              = "${local.name}-alerts"
  kms_master_key_id = aws_kms_key.sns.arn
  tags              = { Name = "${local.name}-alerts" }
}

# Cost monitoring
resource "aws_budgets_budget" "monthly" {
  name              = "${local.name}-monthly-budget"
  budget_type       = "COST"
  limit_amount      = "500"
  limit_unit        = "USD"
  time_unit         = "MONTHLY"
  time_period_start = "2024-01-01_00:00"

  notification {
    comparison_operator        = "GREATER_THAN"
    threshold                  = 80
    threshold_type             = "PERCENTAGE"
    notification_type          = "ACTUAL"
    subscriber_email_addresses = [var.alert_email]
  }

  notification {
    comparison_operator        = "GREATER_THAN"
    threshold                  = 100
    threshold_type             = "PERCENTAGE"
    notification_type          = "FORECASTED"
    subscriber_email_addresses = [var.alert_email]
  }
}

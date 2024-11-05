#!/bin/bash
# job_monitoring_helper.sh
# Purpose: Provide comprehensive job monitoring and log viewing utilities

# Color codes for enhanced readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Job Monitoring and Log Viewing Helper Function
job_monitoring_helper() {
    local job_id="${1:-}"
    
    # Validate job ID
    if [[ -z "$job_id" ]]; then
        echo -e "${RED}Error: Please provide a SLURM job ID${NC}"
        return 1
    }

    # Clear screen and display header
    clear
    echo -e "${BLUE}===== SLURM JOB MONITORING HELPER =====${NC}"
    echo -e "${YELLOW}Job ID: ${job_id}${NC}"
    echo "-------------------------------------------"

    # Section 1: Basic Job Information
    echo -e "${GREEN}1. JOB DETAILS${NC}"
    scontrol show job "$job_id" | grep -E "JobId|JobName|UserId|State|Partition|NodeList"
    echo "-------------------------------------------"

    # Section 2: Real-time Job Status
    echo -e "${GREEN}2. CURRENT JOB STATUS${NC}"
    squeue -j "$job_id"
    echo "-------------------------------------------"

    # Section 3: Resource Utilization
    echo -e "${GREEN}3. RESOURCE UTILIZATION${NC}"
    sstat -j "$job_id" --format=AveCPU,AvePages,AveRSS,MaxRSS
    echo "-------------------------------------------"

    # Section 4: Potential Log Locations
    echo -e "${GREEN}4. POTENTIAL LOG LOCATIONS${NC}"
    echo "SLURM Output Log: slurm-${job_id}.out"
    echo "Custom Logs: ~/logs/[tool]/job_${job_id}/"
    echo "-------------------------------------------"

    # Section 5: Useful Monitoring Commands
    echo -e "${GREEN}5. USEFUL MONITORING COMMANDS${NC}"
    echo -e " ${YELLOW}Live Monitoring:${NC}"
    echo "  watch -n 5 'squeue -j ${job_id}'"
    echo "  seff ${job_id}"
    
    echo -e "\n ${YELLOW}Log Viewing Commands:${NC}"
    echo "  tail -f slurm-${job_id}.out"
    echo "  less slurm-${job_id}.out"
    
    echo -e "\n ${YELLOW}Comprehensive Log Navigation:${NC}"
    echo "  vim ~/logs/[tool]/job_${job_id}/task_*/main_*.log"
    echo "  less ~/logs/[tool]/job_${job_id}/task_*/error_*.log"
    
    echo -e "\n ${YELLOW}Job Cancellation (if needed):${NC}"
    echo "  scancel ${job_id}"
    
    echo -e "\n${BLUE}===== END OF JOB MONITORING HELPER =====${NC}"
}

# Companion Function: Log File Navigator
navigate_job_logs() {
    local job_id="${1:-}"
    local log_base="${HOME}/logs"
    
    if [[ -z "$job_id" ]]; then
        echo -e "${RED}Error: Please provide a SLURM job ID${NC}"
        return 1
    }

    # Find all log directories for the job
    local log_dirs=$(find "${log_base}" -type d -name "job_${job_id}")
    
    if [[ -z "$log_dirs" ]]; then
        echo -e "${RED}No log directories found for job ${job_id}${NC}"
        return 1
    }

    echo -e "${GREEN}Log Directories for Job ${job_id}:${NC}"
    echo "$log_dirs"
    
    echo -e "\n${YELLOW}Suggested Log Exploration Commands:${NC}"
    echo "1. Find all main log files:"
    echo "   find ${log_base} -name \"main_*.log\" -path \"*job_${job_id}*\""
    
    echo -e "\n2. Find all error log files:"
    echo "   find ${log_base} -name \"error_*.log\" -path \"*job_${job_id}*\""
    
    echo -e "\n3. Tail the most recent main log:"
    echo "   tail -n 50 \$(find ${log_base} -name \"main_*.log\" -path \"*job_${job_id}*\" | sort | tail -n 1)"
}

# Usage Examples
: <<'USAGE_EXAMPLES'
# In terminal after job submission
job_monitoring_helper 12345
navigate_job_logs 12345

# Example aliases to add to .bashrc
alias jobmon='job_monitoring_helper'
alias navlogs='navigate_job_logs'
USAGE_EXAMPLES

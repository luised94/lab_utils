#!/bin/bash

# Advanced Logging Functions for Bash Scripts
#
# Script: 002_logging_functions.sh
# Description: A set of functions for consistent logging across Bash scripts
# Author: Your Name
# Date: 2024-10-18

# Function: get_script_name
# Purpose: Extract the full path of the current script
# Parameters: None
# Return: Script path
get_script_name() {
  echo "$(readlink -f "${BASH_SOURCE[0]}")"
}

# Function: get_script_dir
# Purpose: Extract the directory of the current script
# Parameters: None
# Return: Script directory
get_script_dir() {
  echo "$(dirname "$(get_script_name)")"
}

# Function: get_script_basename
# Purpose: Extract the basename of the current script without extension
# Parameters: None
# Return: Script basename
get_script_basename() {
  echo "$(basename "$(get_script_name)" .sh)"
}

# Function: init_logging
# Purpose: Set up logging for a script
# Parameters:
#   $1 - Log file path (optional)
# Return: Log file path
init_logging() {
  local log_file="$1"
  if [[ -z "$log_file" ]]; then
    local script_name="$(get_script_basename)"
    local log_dir="$(get_script_dir)/logs/$(date +%Y-%m)"
    mkdir -p "$log_dir"
    log_file="${log_dir}/$(date +%Y-%m-%d)_${script_name}.log"
  fi
  
  log_system_info "$log_file"
  log_git_info "$log_file"
  
  echo "$log_file"
}

# Function: log_system_info
# Purpose: Log Bash version and system information
# Parameters:
#   $1 - Log file path
# Return: None
log_system_info() {
  local log_file="$1"
  log_message "INFO" "Bash version: ${BASH_VERSION}" "$log_file"
  log_message "INFO" "System: $(uname -a)" "$log_file"
}

# Function: log_git_info
# Purpose: Log git branch and commit hash
# Parameters:
#   $1 - Log file path
# Return: None
log_git_info() {
  local log_file="$1"
  if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    local git_branch=$(git rev-parse --abbrev-ref HEAD)
    local git_hash=$(git rev-parse HEAD)
    log_message "INFO" "Git branch: ${git_branch}" "$log_file"
    log_message "INFO" "Git commit: ${git_hash}" "$log_file"
  else
    log_message "WARNING" "Not in a git repository" "$log_file"
  fi
}

# Function: log_message
# Purpose: Log a message with timestamp and level
# Parameters:
#   $1 - Log level
#   $2 - Message to log
#   $3 - Log file path
# Return: None
log_message() {
  local level="$1"
  local message="$2"
  local log_file="$3"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  local log_entry="[${timestamp}] [${level}] ${message}"
  
  echo "${log_entry}"
  
  if [[ -n "${log_file}" ]]; then
    echo "${log_entry}" >> "${log_file}"
  fi
}

# Function: log_info
# Purpose: Log a message at the INFO level
# Parameters:
#   $1 - Message to log
#   $2 - Log file path
# Return: None
log_info() {
  log_message "INFO" "$1" "$2"
}

# Function: log_warning
# Purpose: Log a message at the WARNING level
# Parameters:
#   $1 - Message to log
#   $2 - Log file path
# Return: None
log_warning() {
  log_message "WARNING" "$1" "$2"
}

# Function: log_error
# Purpose: Log a message at the ERROR level
# Parameters:
#   $1 - Message to log
#   $2 - Log file path
# Return: None
log_error() {
  log_message "ERROR" "$1" "$2"
}

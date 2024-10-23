#!/bin/bash

# File patterns
SLURM_OUTPUT_PATTERN="slurm-*.out"
SLURM_ERROR_PATTERN="slurm-*.err"

# Directory settings
DEFAULT_SLURM_LOG_DIR="slurm_logs"
MAX_LOG_AGE_DAYS=30
MAX_DEPTH_SEARCH=2

# Size limits
MAX_FILE_SIZE_MB=100

# Batch settings
BATCH_SIZE=100

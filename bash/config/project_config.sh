#!/bin/bash
# bash/config/project_config.sh

declare -A PROJECT_CONFIG=(
    [REMOTE_HOST]="luria.mit.edu"
    [REMOTE_USER]="luised94"
    [REMOTE_PATH]="~/data"
    [REQUIRED_DIRS]="documentation fastq"
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"
)

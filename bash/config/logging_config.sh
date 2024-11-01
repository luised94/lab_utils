#!/bin/bash

# bash/config/logging_config.sh
declare -A LOGGING_CONFIG=(
    [RUN_SEPARATOR]="=== New Run ==="
    [TIMESTAMP_FORMAT]="%Y-%m-%d %H:%M:%S"
    [ENTRY_FORMAT]="\n%s (#%d) === %s ===\n"
    [FIRST_RUN_FORMAT]="%s (#1) === %s ==="
    [BUFFER_SIZE]="4096"         # Write buffer size
    [LOCK_BASE_DIR]="/tmp/lab_utils_locks"
    [LOCK_TIMEOUT]=10
    [LOCK_RETRY]=3
    [MAX_MESSAGE_LENGTH]="1024"   # Maximum message length
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"
)


# Critical path protection
declare -rA PROTECTED_PATHS=(
    ["/bin"]=1
    ["/etc"]=1
    ["/dev"]=1
    ["/home"]=1
    ["$HOME"]=1
)

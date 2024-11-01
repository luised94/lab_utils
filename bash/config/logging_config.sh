#!/bin/bash

# bash/config/logging_config.sh
declare -A LOGGING_CONFIG=(
    [RUN_SEPARATOR]="=== New Run ==="
    [TIMESTAMP_FORMAT]="%Y-%m-%d %H:%M:%S"
    [ENTRY_FORMAT]="\n%s (#%d) === %s ===\n"
    [FIRST_RUN_FORMAT]="%s (#1) === %s ==="
    [LOCK_TIMEOUT]="10"           # Seconds to wait for lock
    [LOCK_RETRY]="3"             # Number of retry attempts
    [BUFFER_SIZE]="4096"         # Write buffer size
    [MAX_MESSAGE_LENGTH]="1024"   # Maximum message length
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"
)

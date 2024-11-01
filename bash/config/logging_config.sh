#!/bin/bash

# bash/config/logging_config.sh
declare -A LOGGING_CONFIG=(
    [LOCK_TIMEOUT]="10"           # Seconds to wait for lock
    [LOCK_RETRY]="3"             # Number of retry attempts
    [BUFFER_SIZE]="4096"         # Write buffer size
    [MAX_MESSAGE_LENGTH]="1024"   # Maximum message length
)

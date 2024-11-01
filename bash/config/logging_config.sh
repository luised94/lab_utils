#!/bin/bash

# bash/config/logging_config.sh
declare -A LOGGING_CONFIG=(
    [RUN_SEPARATOR]="=== New Run ==="
    [TIMESTAMP_FORMAT]="%Y-%m-%d %H:%M:%S"
    [ENTRY_FORMAT]="\n%s (#%d) === %s ===\n"
    [FIRST_RUN_FORMAT]="%s (#1) === %s ===\n"
    [BUFFER_SIZE]="4096"         # Write buffer size
    [LOCK_BASE_DIR]="/tmp/lab_utils_locks"
    [LOCK_TIMEOUT]=10
    [LOCK_RETRY]=3
    [MAX_MESSAGE_LENGTH]="1024"   # Maximum message length
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"
)

declare -rA PROTECTED_PATHS=(
    # System directories
    ["^/s?bin(/.*)?$"]=1              # Matches /bin, /sbin and subdirs
    ["^/usr/s?bin(/.*)?$"]=1          # Matches /usr/bin, /usr/sbin and subdirs
    ["^/etc(/.*)?$"]=1                # /etc and all subdirs
    ["^/dev(/.*)?$"]=1                # /dev and all subdirs
    ["^/proc(/.*)?$"]=1               # /proc and all subdirs
    ["^/sys(/.*)?$"]=1                # /sys and all subdirs
    
    # User spaces
    ["^/home$"]=1                     # Just /home
    ["^/root(/.*)?$"]=1               # root's home
    ["^${HOME}$"]=1                   # User's home exact match
    
    # System locations
    ["^/usr/lib(/.*)?$"]=1            # System libraries
    ["^/var/run(/.*)?$"]=1            # Runtime files
    ["^/boot(/.*)?$"]=1               # Boot files
)

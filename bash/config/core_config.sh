# bash/config/core_config.sh
#!/bin/bash

declare -A CORE_CONFIG=(
    # Version Control
    [VERSION]="1.0.0"
    [MIN_BASH_VERSION]="4.2.0"
    # Logging
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_FORMAT]="[%s] [%s] [%s] %s\n"  # timestamp, level, context, message
    # Locking
    [LOCK_TIMEOUT]="30"
    [LOCK_RETRY]="3"
    [LOCK_DIR]="/tmp/lab_utils_locks"
    # Paths
    [PROJECT_ROOT]="$HOME/lab_utils"
    [MODULE_PATH]="$HOME/lab_utils/bash/modules"
    [RUN_SEPARATOR]="=== New Run ==="
    [TIMESTAMP_FORMAT]="%Y-%m-%d %H:%M:%S"
    [ENTRY_FORMAT]="\n%s (#%d) === %s ===\n"
    [FIRST_RUN_FORMAT]="%s (#1) === %s ===\n"
    [BUFFER_SIZE]="4096"         # Write buffer size
    [MAX_MESSAGE_LENGTH]="1024"   # Maximum message length
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

#
#declare -A PATHS=(
#    ["BASE_DIR"]="$HOME/lab_utils"
#    ["FUNCTIONS_DIR"]="bash/functions"
#    ["CONFIG_DIR"]="bash/config"
#    ["SCRIPTS_DIR"]="bash/scripts"
#)
#
#declare -A FILE_PATTERNS=(
#    ["FUNCTIONS"]="*.sh"
#    ["CONFIG"]="*.sh"
#    ["EXCLUDE_PATTERNS"]=(".*" "_*" "test_*")
#)
#
#declare -A LOAD_ORDER=(
#    ["PRIORITY"]=(
#        "logging.sh"
#        "lock_utils.sh"
#        h"
#    )
#    ["OPTIONAL"]=(
#        "experimental.sh"
#        "deprecated.sh"
#    )
#)

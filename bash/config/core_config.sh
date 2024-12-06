#!/bin/bash
# bash/config/core_config.sh

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
    [LOCK_BASE_DIR]="/tmp/lab_utils_locks"
    # Paths
    [PROJECT_ROOT]="$HOME/lab_utils"
    [MODULE_PATH]="$HOME/lab_utils/bash/modules"
    [VERBOSE]="false"
    [RUN_SEPARATOR]="=== New Run ==="
    [TIMESTAMP_FORMAT]="%Y-%m-%d %H:%M:%S"
    [ENTRY_FORMAT]="\n%s (#%d) === %s ===\n"
    [FIRST_RUN_FORMAT]="%s (#1) === %s ===\n"
    [BUFFER_SIZE]="4096"         # Write buffer size
    [MAX_MESSAGE_LENGTH]="1024"   # Maximum message length
)

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
#ÃÄ Step marker
#³  Continuation
#ÀÄ Final step
#û Success
#? Failure

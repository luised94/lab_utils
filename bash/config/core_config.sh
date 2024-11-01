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
)

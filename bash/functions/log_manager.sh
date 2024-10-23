#!/bin/bash

function generate_log_commands() {
    local dir="$1"
    local job_id="$2"
    local time_id="$3"
    
    cat << EOF
View logs:
  Standard output: vim ${dir}/logs/*_${job_id}.out
  Standard error:  vim ${dir}/logs/*_${job_id}.err

Find today's logs:
  find ${dir}/logs -type f -daystart -ctime 0 -name "${LOG_PATTERNS[OUTPUT]}" -exec vim {} +

Find logs by time ID:
  find ${dir}/logs -type f -name "${time_id}*.out" -exec vim {} +

Check file sizes:
  find ${dir}/ -type f -name "${time_id}*" -exec ls -lh {} +
EOF
}

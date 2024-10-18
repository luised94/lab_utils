#!/bin/bash
for file in */*.sh */*/*.sh; do
    if grep -q 'exec >"$OUT_FILE" 2>"$ERR_FILE"' "$file"; then
        sed "$1" '
        /exec >"$OUT_FILE" 2>"$ERR_FILE"/{
            a echo "TASK_START"
        }
        $a echo -e "TASK_END"\n
        ' "$file"
        echo "Modified $file"
    else
        echo "Skipped $file (pattern not found)"
    fi
done
echo "Second for loop starting"
for file in */*.sh */*/*.sh; do
    if grep -q 'OUT_FILE="\${LOG_DIR}/\${timeid}_[^_]*_\${SLURM_ARRAY_JOB_ID}_\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.out' "$file"; then
    sed "$1" '
        s|\(OUT_FILE="\${LOG_DIR}/\${timeid}_[^_]*_\${SLURM_ARRAY_JOB_ID}\)_\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}\(\.out"\)|\1\2|
        s|\(ERR_FILE="\${LOG_DIR}/\${timeid}_[^_]*_\${SLURM_ARRAY_JOB_ID}\)_\${SLURM_JOB_ID}_\${SLURM_ARRAY_TASK_ID}\(\.err"\)|\1\2|
    ' "$file" 
        echo "Modified $file"
    else
        echo "Skipped $file (pattern not found)"
    fi
done

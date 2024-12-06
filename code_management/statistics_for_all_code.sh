#STATUS: REMOVE.
#./code_management/statistics_for_all_code.sh
#!/bin/bash

set -euo pipefail

trap 'echo "Error on line $LINENO. Exit code: $?"; exit 1' ERR

STATS_FILE="repository_stats.txt"
REPO_FILE=$(find . -maxdepth 1 -type f -name "repository_contents.xml")
if [ -f $REPO_FILE ]; then
    echo "File exists, gathering statistics."
else 
    echo "${REPO_FILE} does not exists. Run aggregate_all_code.sh"
    exit 1
fi
echo "Starting repository statistics gathering"


# Gather statistics
echo "Analyzing repository..."
{
    echo "File counts:"
    for file_type in md sh R lua py; do
        count=$(grep -c "<name>.*\.$file_type</name>" "$REPO_FILE")
        echo "$file_type: $count"
    done

    echo "Total lines of code:"
    grep -v '<.*>' "$REPO_FILE" | wc -l

    echo "Top 10 most frequent words:"
    grep -v '<.*>' "$REPO_FILE" | tr '[:upper:]' '[:lower:]' | tr -cs '[:alpha:]' '\n' | sort | uniq -c | sort -rn | head -n 10

    echo "Largest files:"
    grep -n '</file>' "$REPO_FILE" | while IFS=':' read -r line_num _; do
        file_name=$(sed -n "$(($line_num-4))p" "$REPO_FILE" | sed -e 's/.*<name>\(.*\)<\/name>.*/\1/')
        file_size=$(sed -n "$(($line_num-1))p" "$REPO_FILE" | wc -c)
        echo "$file_size $file_name"
    done | sort -rn | head -n 5

    echo "Total final file size:"
    stat -f%z "$REPO_FILE" || echo "Error: Unable to get file size"
} > "$STATS_FILE"

echo "Task completed. Check $REPO_FILE for gathered code and $STATS_FILE for statistics."

# Display stats
cat "$STATS_FILE"

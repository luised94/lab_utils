#!/bin/bash
# bmc_sync_data_folders_to_dropbox.sh
# Purpose: 
#   Sync the data directories with metadata and plots to dropbox.
#   Ensure documentation is shared and available no matter what device I am using.
#   FILES IN THE DROPBOX LOCATION GET OVERWRITTEN!

if [[ -z "${DROPBOX_PATH}" ]]; then
    echo "Error: DROPBOX_PATH not set. Set manually before running."
    echo "DROPBOX_PATH is set via my_config repository."
    echo "Current value: || $DROPBOX_PATH ||"
    exit 1
fi

FROM_DIR="$HOME/data/"
TO_DIR="${DROPBOX_PATH}/Lab/Experiments/ngs"
#excludes=(
#  "fastq"
#  "quality_control"
#  "coverage"
#)
#
## Build exclude flags
#exclude_args=()
#for excl in "${excludes[@]}"; do
#  exclude_args+=(--exclude="$excl/")
#done

if [[ ! -d "${FROM_DIR}" ]]; then
    echo "Error: Experiment directory does not exist: ${EXPERIMENT_DIR}"
    exit 1
fi

if [[ ! -d "${TO_DIR}" ]]; then
    echo "Error: Experiment directory does not exist: ${EXPERIMENT_DIR}"
    exit 1
fi

mapfile -t ALL_DIRECTORIES < <(find "${FROM_DIR}" -mindepth 1 -maxdepth 1 -type d)
declare -a targets

echo "------- Checking directories -----------"
echo "----------------------------------------"
for directory in "${ALL_DIRECTORIES[@]}"; do
  echo "Current directory || ${directory} ||"

  if [[ ! -d ${directory} ]]; then
    echo "Not a directory. Skipping..."
    continue

  fi

  base_directory_name=$(basename "${directory}" )

  if [[ ! "$base_directory_name" =~ ^[0-9]{6}Bel$ ]]; then
    echo "Directory basename not bmc project. Skipping..."
    continue

  fi

  targets+=("$directory")

done
echo "----------------------------------------"

# If no targets, exit cleanly
if [[ ${#targets[@]} -eq 0 ]]; then
  echo "No valid directories found to process."
  exit 0

fi

echo "Found ${#targets[@]} directories to process"

for directory in "${targets[@]}"; do

  echo "$directory"
  base_directory_name=$( basename "${directory}" )
  destination_dir="${TO_DIR}/${base_directory_name}/"
  echo "Basename || ${base_directory_name} ||"
  echo "Destination dir || ${destination_dir} ||"


  echo "Rsync dry run..."
  echo "rsync -nav ${directory}/ $destination_dir"
  rsync -nav \
        --exclude=fastq/ \
        --exclude=quality_control/ \
        --exclude=coverage/ \
        --itemize-changes \
        "${directory}/" "$destination_dir"
  #rsync -avn "${exclude_args[@]}" "$src/" "$dest/"

done
echo "Rsync dry-run complete..."

echo "Execute the rsync command for all of the directory...?"
echo "WARNING: FILES IN THE DROPBOX LOCATION GET OVERWRITTEN!"
read -r -p "Execute all of the above syncs? (y/n): " confirm
echo "Confirm value || $confirm || "
case "$confirm" in
  y|Y)
    echo "Executing all syncs..."
    for directory in "${targets[@]}"; do
      base_directory_name=$( basename "${directory}" )
      destination_dir="${TO_DIR}/${base_directory_name}/"
      rsync -av \
            --exclude=fastq/ \
            --exclude=quality_control/ \
            --exclude=coverage/ \
            --itemize-changes \
            --bwlimit=1000 \
            "${directory}/" "$destination_dir"
      sleep 3
    done
    echo "All syncs completed."
    ;;
  *)
    echo "Aborted. No changes made."
  ;;
esac
echo "Finished rsync procedure..."

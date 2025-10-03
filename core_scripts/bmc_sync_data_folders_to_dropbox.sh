#!/bin/bash

if [[ -z "${DROPBOX_PATH}" ]]; then
    echo "Error: DROPBOX_PATH not set. Set manually before running."
    echo "DROPBOX_PATH is set via my_config repository."
    echo "Current value: || $DROPBOX_PATH ||"
    exit 1
fi

FROM_DIR="$HOME/data/"
TO_DIR="${DROPBOX_PATH}/Lab/Experiments/ngs"

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

for directory in "${ALL_DIRECTORIES[@]}"; do
  echo "----------------------------------------"

  echo "Current directory || ${directory} ||"

  if [[ ! -d ${directory} ]]; then
    echo "Not a directory. Skipping..."
    echo "----------------------------------------"
    continue
  fi

  base_directory_name=$(basename "${directory}" )

  if [[ ! "$base_directory_name" =~ ^[0-9]{6}Bel$ ]]; then
    echo "Directory basename not bmc project. Skipping..."
    echo "----------------------------------------"
    continue
  fi

  destination_dir="${TO_DIR}/${base_directory_name}/"
  echo "Basename || ${base_directory_name} ||"
  echo "Destination dir || ${destination_dir} ||"


  echo "Syncing directory..."
  echo "rsync -nav ${directory}/ $destination_dir"
  rsync -nav \
        --exclude=fastq/ \
        --exclude=quality_control/ \
        --exclude=coverage/ \
        --itemize-changes \
        "${directory}/" "$destination_dir"

  targets+=("$directory")
  echo "----------------------------------------"
done
#for directory in "${ALL_DIRECTORIES[@]}"; do

# If no targets, exit cleanly
if [[ ${#targets[@]} -eq 0 ]]; then
  echo "No valid directories found to process."
  exit 0
fi

#echo "Rsync dry-run complete..."
#
#
#echo
#echo "Execute the rsync command for all of the directory...?"
#read -r -p "Execute all of the above syncs? (y/n): " confirm
#echo "Confirm value || $confirm || "
#case "$confirm" in
#  y|Y)
#    echo "Executing all syncs..."
#    for directory in "${ALL_DIRECTORIES[@]}"; do
#      basename="$(basename "$src")"
#      dest="$destination_base/$basename"
#      echo "Syncing: $basename"
#      rsync -av "${exclude_args[@]}" "$src/" "$dest/"
#    done
#    echo "All syncs completed."
#    ;;
#  *)
#    echo "Aborted. No changes made."
#  ;;
#esac
#echo "Potato"


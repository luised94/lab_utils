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

for directories in "${ALL_DIRECTORIES[@]}"; do
  echo "----------------------------------------"

  base_directory_name=$(basename "${directories}" )
  echo "Current directory|| ${directories} ||"
  echo "Basename|| ${base_directory_name} ||"

  if [[ ! "$base_directory_name" =~ ^[0-9]{6}Bel$ ]]; then
    echo "Directory basename not bmc project. Skipping..."
    continue
  fi

  echo "Syncing directory..."
  echo "rsync -nav ${directories}/ ${TO_DIR}/${base_directory_name}/"
  rsync -nav "${directories}/" --exclude=fastq --exclude=quality_control --exclude=coverage "${TO_DIR}/${base_directory_name}/"

  echo "----------------------------------------"
done

#!/bin/bash

DIRECTORY="$HOME/data/250930Bel/fastq"
NUMBER_OF_LANES=3
IS_PAIRED_END=1
FASTQ_FILES_PER_LANE=
# Example of filename: 250930Bel_D25-12496-2_2_sequence.fastq
# Check the fastq files using just _ and -_?
IDX_FOR_ID_PART=2
IDX_FOR_LANE_PART=3
IDX_FOR_READ_PAIR_PART=4

if [[ IS_PAIRED_END -eq 1 ]];then
  FASTQ_FILES_PER_LANE=2
  EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE=$(( NUMBER_OF_LANES * FASTQ_FILES_PER_LANE))

else
  FASTQ_FILES_PER_LANE=1
  EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE=$(( NUMBER_OF_LANES * FASTQ_FILES_PER_LANE))

fi


echo "Expected number of fastq files: $EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE"
mapfile -t fastq_files < <(find "$DIRECTORY" -type f -name "*.fastq")

echo "Getting unique_ids"
readarray -t unique_ids < <(
  for fastq_file in "${fastq_files[@]}"; do
    if [[ "$fastq_file" =~ unmapped ]]; then
      echo "Skipping unmapped fastq file"
      continue
    fi
    # Split filename components
    IFS='_-' read -ra parts <<< "${fastq_file##*/}"
    # Extract ID from standardized position (3rd component)
    printf "%s\n" "${parts[$IDX_FOR_ID_PART]}"
  done | sort -u | grep -v "unmapped"
)

echo "Getting unique lanes..."
readarray -t unique_lanes < <(
  for fastq_file in "${fastq_files[@]}"; do
    if [[ "$fastq_file" =~ unmapped ]]; then
      echo "Skipping unmapped fastq file"
      continue
    fi
    # Split filename components
    IFS='_-' read -ra parts <<< "${fastq_file##*/}"
    # Extract ID from standardized position (3rd component)
    printf "%s\n" "${parts[$IDX_FOR_LANE_PART]}"
  done | sort -u | grep -v "unmapped"
)

echo "Getting unique read pairs..."
readarray -t unique_read_pairs_id < <(
  for fastq_file in "${fastq_files[@]}"; do
    if [[ "$fastq_file" =~ unmapped ]]; then
      echo "Skipping unmapped fastq file"
      continue
    fi
    # Split filename components
    IFS='_-' read -ra parts <<< "${fastq_file##*/}"
    # Extract ID from standardized position (3rd component)
    printf "%s\n" "${parts[$IDX_FOR_READ_PAIR_PART]}"
  done | sort -u | grep -v "unmapped"
)

printf '%s\n' "${unique_ids[@]}" | xargs -n6 | sed 's/^/  /' | column -t
printf '%s\n' "${unique_lanes[@]}" | xargs -n6 | sed 's/^/  /' | column -t
printf '%s\n' "${unique_read_pairs_id[@]}" | xargs -n6 | sed 's/^/  /' | column -t

for id in "${unique_ids[@]}"; do
  echo "--------------------"
  echo "Processing ID: $id"
  # Find files using array and glob pattern
  #files=( *"${id}"*_sequence.fastq )
  mapfile -t files < <( find "$DIRECTORY" -type f -name "*${id}*" )

  if [[ ${#files_with_id[@]} -eq 0 ]]; then
    echo "ERROR: No files found for ID: $id"
    continue
  fi

  if [[ ! ${#files_with_id[@]} -eq $EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE ]]; then
    echo -e "[ERROR]: Found ${#files_with_id[@]} for $id\nExpected $EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE"
    continue
  fi
  echo " ID '$id' has expected number of fastq files."

  for lane in "${unique_lanes[@]}"; do
    mapfile -t files_in_lane < <( find "$DIRECTORY" -type f -name "*${id}-${unique_lanes}*" )
  done


done

#echo "${unique_ids[@]}"
#echo "${fastq_file[@]}" | awk -F/ '{print $NF}'
#cat ${fastq_files[6]} | head -n 8
#cat ${fastq_files[7]} | head -n 8
#cat ${fastq_files[6]} | wc -l
#cat ${fastq_files[7]} | wc -l
#for fastq_file in "${fastq_files[@]}"; do
#  if [[ "$fastq_file" =~ unmapped ]]; then
#    echo "Skipping unmapped fastq file"
#    continue
#  fi
#  echo "$fastq_file"
#
#done

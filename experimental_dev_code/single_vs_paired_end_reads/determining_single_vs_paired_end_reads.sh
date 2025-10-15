#!/bin/bash

DIRECTORY="$HOME/data/250930Bel/fastq"
NUMBER_OF_LANES=3
IS_PAIRED_END=1

if [[ IS_PAIRED_END -eq 1 ]];then
  EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE=$(( NUMBER_OF_LANES * 2))
else
  EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE=$(( NUMBER_OF_LANES * 1))
fi
echo "Expected number of fastq files: $EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE"

mapfile -t fastq_files < <(find "$DIRECTORY" -type f -name "*.fastq")

cat ${fastq_files[6]} | head -n 8
cat ${fastq_files[7]} | head -n 8

cat ${fastq_files[6]} | wc -l
cat ${fastq_files[7]} | wc -l

# Example of filename: 250930Bel_D25-12496-2_2_sequence.fastq
for fastq_file in "${fastq_files[@]}"; do
  if [[ "$fastq_file" =~ unmapped ]]; then
    echo "Skipping unmapped fastq file"
    continue
  fi
  echo "$fastq_file"

done

readarray -t unique_ids < <(
  for fastq_file in "${fastq_files[@]}"; do
    if [[ "$fastq_file" =~ unmapped ]]; then
      echo "Skipping unmapped fastq file"
      continue
    fi
    # Split filename components
    IFS='_-' read -ra parts <<< "${fastq_file##*/}"
    # Extract ID from standardized position (3rd component)
    printf "%s\n" "${parts[2]}"
  done | sort -u | grep -v "unmapped"
)

printf '%s\n' "${unique_ids[@]}" | xargs -n6 | sed 's/^/  /' | column -t
for id in "${unique_ids[@]}"; do
  if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: No files found for ID: $id"
    continue
  fi

  if [ ! ${#files[@]} -eq 2 ]; then
    echo -e "[ERROR]: Found ${#files[@]} for $id\nExpected 2"
    continue
  fi

done
#echo "${unique_ids[@]}"
#echo "${fastq_file[@]}" | awk -F/ '{print $NF}'

for id in "${unique_ids[@]}"; do
  echo "--------------------"
  echo "Processing ID: $id"
  # Find files using array and glob pattern
  #files=( *"${id}"*_sequence.fastq )
  mapfile -t files < <( find "$DIRECTORY" -type f -name "*${id}*" )

  if [[ ${#files[@]} -eq 0 ]]; then
    echo "ERROR: No files found for ID: $id"
    continue
  fi

  if [[ ! ${#files[@]} -eq $EXPECTED_NUMBER_OF_FASTQ_FILES_PER_SAMPLE ]]; then
    echo -e "[ERROR]: Found ${#files[@]} for $id\nExpected 2"
    continue
  fi
  echo "Files have expected number of fastq files."
done

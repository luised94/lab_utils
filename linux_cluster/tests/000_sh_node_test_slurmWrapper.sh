#!/bin/bash


if [ $# -ne 3 ]; then
	echo "Usage: $0 <array_number> <directory_to_process> <script_name>"
	echo "Array number is an integer depending on the number of array tasks to create (--array="
	echo 'Directory is directory to process with / ( eg 240304Bel/ )'
	echo 'script_name is absolute path to some lab_utils script ( eg ~/lab_utils/next_generation_sequencing/000/000.sh )'
	exit 1
fi

echo "sbatch --array=$1 $2 $3"

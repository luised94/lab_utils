#!/bin/bash

accessions=("GCA_000146045.2" "GCA_000005845.2" "GCA_000001405.29" "GCA_002163515.1")
download_dirs=("${accessions[@]/%/_genome}")

for download_dir in "${download_dirs[@]}"; do
	echo ${download_dir}
done

#!/bin/bash


find . -type f -name "*_refgenome.fna" -exec sh -c 'bowtie2-build "$0" "${0%_refgenome.fna}_index"' {} \;
#find . -type f -name "*_refgenome.fna" -exec sh -c 'echo "$0" "$0%_refgenome.fna}_index"' {} \;

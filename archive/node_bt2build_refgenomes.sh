#!/bin/bash
#USAGE: From ~/data/REFGENS/, run './node_ bt2build_refgenomes.sh'
#
refgenome_dir=~/data/REFGENS
exec > "${refgenome_dir}/logs/$(date '+%Y-%m-%d-%M-%S')_indexing.out" 2>&1
echo "START TIME: $(date '+%Y-%m-%d-%M-%S')"
echo "From dir ${refgenome_dir}"

module purge
module load gnu/5.4.0
module load bowtie2/2.3.5.1

mapfile -t genome_paths < <(find "$refgenome_dir" -type f -name "*_refgenome.fna")

#genome_path=${genome_paths[$SLURM_ARRAY_TASK_ID-1]}
for genome_path in "${genome_paths[@]}"; do
	echo "Starting indexing for $genome_path"
	echo "$genome_path" "${genome_path%_refgenome.fna}_index"
	bowtie2-build "$genome_path" "${genome_path%_refgenome.fna}_index"
	echo "Indexing completed"
done 

echo "END TIME: $(date "+%Y-%m-%d-%M-%S")"


# config.sh.template
# Copy to config.sh and customize. config.sh is gitignored.

#==============================
# Reference Files
#==============================
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="${GENOME_DIR}/SaccharomycescerevisiaeS288C_index"
BLACKLIST_BED_FILE="$HOME/data/feature_files/20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed"

#==============================
# fastp
#==============================
MAXIMUM_N_BASE_COUNT=5
MAXIMUM_UNQUALIFIED_BASE_PERCENT=50
COMPLEXITY_WINDOW_SIZE=4
OVERREPRESENTATION_SAMPLING=50
MINIMUM_BASE_QUALITY=20
MINIMUM_READ_LENGTH=25
FASTP_OUT_PREFIX="fastpfiltered"

#==============================
# bowtie2
#==============================
ALIGNMENTS_TO_REPORT=1
EXTENSIONS_TO_TRY=15
SETS_OF_SEEDS=2
MAX_PENALTY_MISMATCH=4
NON_NUCLEOTIDE_PENALTY=1
BAM_OUT_SUFFIX="S288C_sorted"
BLFILTERED_OUT_SUFFIX="blFiltered"

#==============================
# bamCoverage
#==============================
BIN_SIZE=10
EFFECTIVE_GENOME_SIZE=12157105
MINIMUM_MAPPING_QUALITY=20
NORM_METHOD="RAW"
VALID_NORM_METHODS=("RAW" "RPKM" "CPM" "BPM" "RPGC")

#==============================
# Naming Conventions
#==============================
FASTQ_FILEPATH_PATTERN="consolidated*.fastq"
MANIFEST_FILENAME="consolidated_reads_manifest.tsv"

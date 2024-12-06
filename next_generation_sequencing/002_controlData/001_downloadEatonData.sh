#STATUS:
#!/bin/bash
#DESCRIPTION: Download data from Eaton 2010 paper as control and comparison. 
#USAGE: 001_downloadEatonData.sh <dir>
#NOTE: Only meant to be run once to download control data from Eaton 2010 paper.
#TODO: Provide arguments for files to concatenate and accession to generalize the script. 

#BIOPROJECT_ACCESSION=PRJNA117641
#Eaton2010 
#Title ORC precisely positios nucleosomes at origins of replication
#DOWNLOAD ORC WT samples in G2 arrest.
DOWNLOAD_TO_DIR="$1"
OUTPUT_DIR="${HOME}/data/${DOWNLOAD_TO_DIR}"
mkdir -p $OUTPUT_DIR
BASE_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"

mapfile -t FILES_TO_CONCATENATE < <(printf "%s\n" WT-G2-ORC-rep{1,2}.fastq.gz)
mapfile -t ACCESSION < <(printf "%s\n" SRR0344{75,76})

FILE_INDEX=0
for FILE in ${FILES_TO_CONCATENATE[@]}; do 
    OUTPUT_FILE=${OUTPUT_DIR}$FILE
    URL_TO_DOWNLOAD=${BASE_URL}${ACCESSION[FILE_INDEX]:0:6}/${ACCESSION[FILE_INDEX]}/${ACCESSION[FILE_INDEX]}.fastq.gz 
    printf "%s | %s | %s\n" "$FILE_INDEX" "$OUTPUT_FILE" "$URL_TO_DOWNLOAD"
    curl -I $URL_TO_DOWNLOAD 
  #  wget --output-document=$OUTPUT_FILE $URL_TO_DOWNLOAD
  #  cat $OUTPUT_FILE >> $OUTPUT_DIR/nnNnH.fastq.gz
    ((FILE_INDEX++))
done 

gunzip ${OUTPUT_DIR}nnNnH.fastq.gz
#echo "gunzip ${OUTPUT_DIR}nnNnH.fastq.gz"

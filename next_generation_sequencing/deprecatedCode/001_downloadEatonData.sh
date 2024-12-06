#STATUS:
#!/bin/bash
#BIOPROJECT_ACCESSION=PRJNA117641
#Eaton2010 
#Title ORC precisely positios nucleosomes at origins of replication
#DOWNLOAD ORC WT samples in G2 arrest.
DOWNLOAD_TO_DIR="$1"
OUTPUT_DIR="${HOME}/data/${DOWNLOAD_TO_DIR}"

FILENAMES=(WT-G2-ORC-rep1.fastq.gz WT-G2-ORC-rep2.fastq.gz)
WEBSITES=(ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034475/SRR034475.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034476/SRR034476.fastq.gz)

FILE_INDEX=0
for FILE in ${FILENAMES[@]}; do 
  OUTPUT_FILE=${OUTPUT_DIR}/$FILE
  printf "%s | %d | %s\n" "$FILE_INDEX" "$OUTPUT_FILE" "${WEBSITES[FILE_INDEX]}"
  wget --output-document=$OUTPUT_FILE ${WEBSITES[COUNTER]}
  cat $OUTPUT_FILE >> $OUTPUT_DIR/nnNnH.fastq.gz
  rm $OUTPUT_FILE
  ((FILE_INDEX++))
done 

gunzip $OUTPUT_DIR/nnNnH.fastq.gz
#echo "gunzip $OUTPUT_DIR/nnNnH.fastq.gz"

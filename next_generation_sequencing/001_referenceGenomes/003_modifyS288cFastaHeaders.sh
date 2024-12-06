#STATUS: REMOVE.
# Description: Renames the fasta hearders for the S288c genome to fit USCS specification.
# Usage: I ran the lines individually on the command line as a way to test. 
# Not idempotent so cant continually test. Run only once.
#Find the S288C genome to modify the headers to UCSC standards
#TODO: Determine if the while loop works as it looks more efficient.
REFGENOME_DIR="$HOME/data/REFGENS"
PATH_TO_S288C_GENOME=$(find $REFGENOME_DIR -type f -name "*S288C_refgenome.fna")
#Create a backup file. 
BACKUP_FILE="${PATH_TO_S288C_GENOME%_refgenome.fna}_backup.fna"
cat $PATH_TO_S288C_GENOME > $BACKUP_FILE

# Process S288C genome file: Convert chromosome names, reformat headers, and remove extra information
awk '/^>/ {gsub(/chromosome/, "chr", $6); printf(">%s%s\n", $6, $7)} !/^>/ {print $0}' "$BACKUP_FILE" | sed 's/,.*//' > $PATH_TO_S288C_GENOME

#while IFS= read -r line; do
#  if [[ $line == ">"* ]]; then
#    echo "${line/chromosome/chr}" | cut -d',' -f1
#  else
#    echo "$line"
#  fi
#done < "$BACKUP_FILE" > "$PATH_TO_S288C_GENOME"
#To recreate the original from backup 
#cat $BACKUP_FILE > $PATH_TO_S288C_GENOME

#Find the S288C genome to modify the headers to UCSC standards
REFGENOME_DIR="$HOME/data/REFGENS"
PATH_TO_S288C_GENOME=$(find $REFGENOME_DIR -type f -name "*S288C_refgenome.fna")
#Create a backup file. 
BACKUP_FILE="${PATH_TO_S288C_GENOME%_refgenome.fna}_backup.fna"
cat $PATH_TO_S288C_GENOME > $BACKUP_FILE

awk '/^>/ {gsub(/chromosome/, "chr", $6); printf(">%s%s\n", $6, $7)} !/^>/ {print $0}' "$BACKUP_FILE" | sed 's/,.*//' > $PATH_TO_S288C_GENOME

#To recreate the original from backup 
#cat $BACKUP_FILE > $PATH_TO_S288C_GENOME

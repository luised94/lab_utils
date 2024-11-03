# 

#

# Steps
Run srun according to email from bmc. 
```bash
srun rsync -av /net/bmc-pub17/data/bmc/public/Bell/241007Bel ~/data/241007Bel/
```

find . -type f -name "*.fastq" -exec mv {} . \; || {
        echo "Failed to move FASTQ files"
        cd "$current_dir"
        return 1
}

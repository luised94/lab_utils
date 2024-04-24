# Running rsync through R ----
#The following requires WSL to be installed on Windows computer 
#
# dir_to_rsync_from  <- "../ngs-pipeline/scripts/R-files/*"
# rsync_dest <- "./scripts/"
# rsync_command_dry <- 'wsl rsync --stats -nv'
# rsync_command <- 'wsl rsync --stats -v'
#
# #Check the string
# paste(rsync_command_dry, dir_to_rsync_from, rsync_dest)
#
# #Check the command output
# system(paste(rsync_command_dry, dir_to_rsync_from, rsync_dest))

#Transfer if it looks good
# system(paste(rsync_command, dir_to_rsync_from, rsync_dest))

# #Transfer files to luria server 
# rsync_command_dry <- 'wsl rsync --stats -nvr'
# rsync_command <- 'wsl rsync --stats -vr'
# dir_to_rsync_from  <- stringr::str_replace(paste0(getwd(), "/scripts/"), pattern = "C:", replacement = "/mnt/c")
# rsync_dest <- "luised94@luria.mit.edu:/home/luised94/data/R-scripts/"
# #For this you can run from command line. It requires password. Just paste output from R console to the terminel
# paste(rsync_command_dry, dir_to_rsync_from, rsync_dest)
# # system(paste(rsync_command_dry, dir_to_rsync_from, rsync_dest))

#Run rsync on terminal ----
#Sync to dropbox
wsl rsync --stats -nrv /mnt/c/Users/Luis/Projects/* '/mnt/c/Users/Luis/Dropbox (MIT)/Lab/Code/Projects/'
#Rsync my Zotero folder
wsl rsync -nrv /mnt/c/Users/Luis/Zotero/* '/mnt/c/Users/Luis/Dropbox (MIT)/Zotero/'
wsl rsync -nrv --update /mnt/c/Users/Luis/Zotero/* '/mnt/c/Users/Luis/Dropbox (MIT)/Zotero/'

#Inverted the order to sync from the cluster to the local machine
wsl rsync --stats -nvr --update luised94@luria.mit.edu:/home/luised94/data/R-scripts/* /mnt/c/Users/Luis/Projects/working-on-a-cluster/scripts/cluster-modified/

#Transfering and updating to luria
wsl rsync --stats -nv --update /mnt/c/Users/Luis/Projects/working-on-a-cluster/scripts/* luised94@luria.mit.edu:/home/luised94/data/rscripts/

#Count the lines that have fastq in it. 
srun rsync -nav /net/bmc-pub17/data/bmc/public/Bell/221024Bel /home/luised94/data/221024Bel_CHIP/data/fastq-files | grep -e ".fastq"

find -type d -name "*fastq*" | grep "/f"


#Run for files after rsync ----
find -type f -exec dos2unix -k -o {} \;
find -type f -exec chmod +x {} \;
#----
#Downloading Eaton data 
wget --output-document=./fastq-files/WT-G2-ORC-rep1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034475/SRR034475.fastq.gz
wget --output-document=./fastq-files/WT-G2-ORC-rep2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034476/SRR034476.fastq.gz 

#Download BMC data 
srun rsync -nav /net/bmc-pub17/data/bmc/public/Bell/221024Bel /home/luised94/data/221024Bel_CHIP/data/fastq-files

#Move the folder content up one directory
mv 221024Bel/* . 
rmdir 221024Bel


sacct -A luised94 --format=JobName,Account,AllocNodes,AllocCPUs,AllocTres,AveDiskRead,Elapsed,JobID,ReqMem | awk '$1~/sbatch/'

test_var="$(echo date +%g_%m_%d_%H_%M_%S)"
${test_var}

echo Execution Time $(${test_var}) > test.out

sed -i '2 i '"$(echo Execution Time $($test_var))"'' test.out


myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"
#Moving the raw fastq data from the Dropbox folder to the orc4r_suppressor_screen folder and processing certain files to label the files with consistent naming. 


#To get target directory, use $echo "$(pwd)" > working_dict.txt; then copy and paste as target_directory variable

target_directory="/mnt/c/Users/liusm/Documents/Test_Directory/orc4r_suppressor_screen/data/seqs"

#From the Dropbox directory C:\Users\liusm\Dropbox (MIT)\Sequencing Data for ORC4 Suppressor Screen
find -name "*sequence.fastq" | xargs -n 1 cp $target_directory 


#Processing fastaseq_id file so that it is tabular data and can be used by cut and column. 
#Use sed s/pattern/replacement to serially replace text in the file to turn it into a table. 
#fasta_id.txt can be used with cut and column.

cat fastafile_id.txt | sed s/":"/""/g | sed s/"  "/" "/g | sed s/" "/"\t"/g > fasta_id.txt


#Use awk twice on fasta_id.txt to add column with integer to identify number of file as well the connection between the files. Output first nine samples then append the last nine samples. S1 and C1 should be compared at after variant calling, S2 and C2 should be compared, etc, etc. 

awk 'BEGIN{ s = 0}; {s += 1}; s<= 9 {print s "\t" "S"s "\t" $1 "\t" $2 "\t" $3 }' fasta_id.txt > fastq_id.txt

awk 'BEGIN{ s = 0}; {s += 1}; s > 9 {print s "\t" "C"s-9 "\t" $1 "\t" $2 "\t" $3 }' fasta_id.txt >> fastq_id.txt

#Use the generated fastq_id.txt file to label the raw fastq data. This file is a few directories up. Extract the second column of fastq_id.txt that has the corresponding sample labeling.   

sample_id=($(cut -f 2 fastq_id.txt))
echo ${sample_id[@]}

#Go to directory with the raw sequencing.
cd ./data/seqs/

#Output file with name of raw fastq files and extract to create array.
find -name "*sequence.fastq" > raw_fastq.txt

sample_file=($(cut -f 1 raw_fastq.txt))
echo ${sample_file[@]}

#Use for loop to iterate through the array with the name of the fastq files. Create counter to go through the sample_id array. Use mv and awk to replace 191209Bel with the sample label. 

COUNTER=0
for file in ${sample_file[@]}; do
   mv -i -T $file "${sample_id[COUNTER]}""$(awk '{ sub(/.*191209Bel/, ""); print }' <<< $file)"
   let COUNTER++;
done

#Other example 
for file in *.jpeg; do
    mv -- "$file" "$(basename $file .jpeg).jpg"
done


#!/bin/bash
cd /home/luised94/data/221024Bel_CHIP/data/fastq-files/
pwd

outputfiles=(WT-G2-ORC-rep1.fastq.gz WT-G2-ORC-rep2.fastq.gz)
websites=(ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034475/SRR034475.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034476/SRR034476.fastq.gz)

COUNTER=0
for file in ${outputfiles[@]};
do 
  echo $file
  echo $COUNTER 
  echo ${websites[COUNTER]}
  wget --output-document=$file ${websites[COUNTER]}
  cat $file >> nnNnH.fastq.gz
  rm $file
  let COUNTER++
done 

# #Check that file is right size
# ls -l --block-size=M
#Do this after files have been renamed
gunzip nnNnH.fastq.gz
#cat 221024Bel_CHIP.txt | sed s/"    "/""/g | sed s/":"/""/g | sed s/"  "/" "/g | sed s/" "/"\t"/g 
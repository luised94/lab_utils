# README for Example Workflow of CHIP-seq Data

> Contains information about the files used for the example analysis, commands used and the main resource. 

>Take a look at ./preliminary-steps-chip-seq-analysis.md to find out what you should install first 

### Create virtual environment with python 2 and create directories

'$virtualenv -p 'location/of/python2.7' 'name-of-directory'

'$mkdir genome-files fastq-files fastqc-files sam-files peak-files'


### Download and index reference genome files from NCBI

> Download the files from NCBI directory. Depending on organism and strain, the download may change. Consult scientific literature to make sure you have downloaded the correct one. 

> If you are able to use and validate this file, you can use it for other experiment as long as it is applicable. Just make sure to note when it was downloaded and how.

'$wget --directory-prefix=genome-files/R64/ -e robots=off -nd -nH --reject="index.html*" --mirror --cut-dirs=2 --recursive --no-parent --convert-links https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/146/045/GCA_000146045.2_R64/ '

'$ wget --directory-prefix=genome-files/w303-asm/ -e robots=off -nd -nH --reject="index.html*" --mirror --cut-dirs=2 --no-parent --convert-links https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/163/515/GCA_002163515.1_ASM216351v1/ '

> Can use command provided by NCBI but the output is different

'$curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/GCF_000146045.2/download?include_annotation_type=GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA&filename=GCF_000146045.2.zip" -H "Accept: application/zip"'

> Can use the reference genome provided by Bioconductor 

'> install.packages("BiocManager")'
'> BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")'
'> genome <- BSgenome.Scerevisiae.UCSC.sacCer3'


> Download reference files from 'Rossi, M.J., Kuntala, P.K., Lai, W.K.M. et al. A high-resolution protein architecture of the budding yeast genome. Nature 592, 309–314 (2021). https://doi.org/10.1038/s41586-021-03314-8 ' Will probably have to adjust labels compared to reference genome downloaded. 

'$svn export https://github.com/CEGRcode/2021-Rossi_Nature.git/trunk/02_References_and_Features_Files '

> You can download the SGD_features directly from yeastmine.org. A python file can be generated to repeat the script programmatically.

> Download Hawkins paper origin information to same directory

'$wget -O "hawkins-origins-timing.xlsx" https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx ' 


> Rename files to refer to them easily. This is optional and you should remember it when running future commands. In this case, I used find and rename to replace some of the filenames so that they do not overlap but still contain some useful information. Use -n option to display the files and there replacement name to tune the commands. 

'$find  -iname "GCA*R64_rna*" | rename -n 's/GCA.*_/R64_rna_/' '
'$find  -iname "GCA*R64_cds*" | rename -n 's/GCA.*_/R64_cds_/' '
'$find  -iname "GCA*R64_feature*" | rename -n 's/GCA.*_/R64_feature_/' '
'$find  -iname "GCA*R64*" | rename 's/GCA.*_/R64_/' '

'$find  -iname "GCA*ASM*" | rename 's/GCA.*_/w303-asm_/' '

'$md5sum path/to/file/R64_genomic.fna'

> Convert excel files to csv to load into R later on. 
> Unzip and index the file. 

'$gzip -d ./genome-files/R64/R64_genomic.fna.gz ''

'$bowtie-build ./genome-files/R64/R64_genomic.fna ./genome-files/R64/R64_genomic_index '

'$bowtie-build ./genome-files/w303-asm/w303-asm_genomic.fna ./genome-files/w303-asm/w303-asm_genomic_index ' 

### Download and QC read files 

> Depending on the analysis that you are doing, you may have to download the files from a core or institutional server. 

> Here we are using a file from an already published paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2854390/). Find the experiment accession (SRX****) from the GEO website (https://www.ncbi.nlm.nih.gov/geo/) and find the fastq data from the ENA website (https://www.ebi.ac.uk/ena/browser/home)

'$wget --output-document=./fastq-files/WT-G2-ORC-rep1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034475/SRR034475.fastq.gz '

'$wget --output-document=./fastq-files/WT-G2-ORC-rep2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034476/SRR034476.fastq.gz '

> Unzip the files 

'$find -iname "WT-G2-ORC*" | xargs -n 1 gzip -d '

> Run fastqc program and assess the quality of the reads. When aligning, you usually want to map only the high-quality bases. To better understand the fastqc program, go to http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010 and to https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

'$find -iname "WT-G2-ORC*" | xargs -n 1 fastqc --outdir ./fastqc-files/ '

> Map the fastq file to the reference genome. View the out file to check important statistics. Use the adapter sequences found in the fastqc data to normalize the amount of reads mapped.  

'$bowtie ./genome-files/w303-asm/w303-asm_genomic_index -q ./fastq-files/WT-G2-ORC-rep1.fastq -v 1 -m 1 -3 9 -S 2> ./sam-files/WT-G2-ORC-rep1.out > ./sam-files/WT-G2-ORC-rep1.sam '

'$bowtie ./genome-files/w303-asm/w303-asm_genomic_index -q ./fastq-files/WT-G2-ORC-rep2.fastq -v 1 -m 1 -3 9 -S 2> ./sam-files/WT-G2-ORC-rep2.out > ./sam-files/WT-G2-ORC-rep2.sam '

> Use macs to call peaks in each file. You dont need an input control file but it does help.

macs -t ./sam-files/WT-G2-ORC-rep1.sam --format SAM --gsize 10000000 --name "WT-G2-ORC-rep1" - --outdir='./peak-files/' --bw 400 --keep-dup 1 --bdg --single-profile --diag &> ./peak-files/WT-G2-ORC-rep1-MACS.out 

macs -t ./sam-files/WT-G2-ORC-rep2.sam --format SAM --gsize 10000000 --name "WT-G2-ORC-rep2" - --outdir='./peak-files/' --bw 400 --keep-dup 1 --bdg --single-profile --diag &> ./peak-files/WT-G2-ORC-rep2-MACS.out 



>Link to resource: https://www.biologie.ens.fr/~mthomas/other/chip-seq-training/#


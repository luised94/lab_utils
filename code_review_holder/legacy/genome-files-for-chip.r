library(BSgenome)
library(GenomeInfoDb)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(GenomicRanges)
library(GenomicAlignments)
library(ComplexHeatmap)
library(circlize)
setwd(file.path("C:", "Users", "Luis", "virtual-environments", "chip-training-yeast-py2", "genome-files", "R64"))



genomefilePath = paste(getwd(), "/R64_genomic.fna", sep="")

#Load genome file from folder and from BSgenome package
sacCer <- readDNAStringSet(genomefilePath)
sacCer1 <-BSgenome.Scerevisiae.UCSC.sacCer3 

#At this point seqlevelStyle(sacCer) doesnt work
seqlevelsStyle(sacCer)
seqlevelsStyle(sacCer1)

#Can use BSgenome to set names to UCSC scheme
names(sacCer1)
names(sacCer) <- names(sacCer1)[1:16]

#Check that genome still works
sacCer[names(sacCer)[1]]
sacCer[names(sacCer)[7]]
seqlevelsStyle(sacCer)

#Creating tiled chromosome as a GRange object 

cere_chrs = getChromInfoFromUCSC('sacCer3')

cere_chrs = subset(cere_chrs, grepl('chrVII$', chrom))

seqlengths = with(cere_chrs, setNames(size, chrom))

tilling_window = tileGenome(seqlengths, tilewidth=1000)

tilling_window = unlist(tilling_window)
tilling_window

seq = getSeq(sacCer1, tilling_window)

# calculates the frequency of all possible dimers 
# in our sequence set
nuc = oligonucleotideFrequency(seq, width = 2)

# converts the matrix into a data.frame
nuc = as.data.frame(nuc)

# calculates the percentages, and rounds the number
nuc = round(nuc/1000,3)

# counts the number of reads per tilling window 
# for each experiment
so = summarizeOverlaps(tilling_window, bam_files)

# converts the raw counts to cpm values
counts  = assays(so)[[1]]
cpm     = t(t(counts)*(1000000/colSums(counts)))

# because the cpm scale has a large dynamic range
# we transform it using the log function
cpm_log = log10(cpm+1)

gc = cbind(data.frame(cpm_log), GC = nuc['GC'])


library(Rsamtools)
library(GenomicRanges)

# read CpGi data set
filePath=system.file("extdata",
                     "cpgi.hg19.chr21.bed",
                     package="compGenomRData")
cpgi.df = read.table(filePath, header = FALSE,
                     stringsAsFactors=FALSE) 
# remove chr names with "_"
cpgi.df =cpgi.df [grep("_",cpgi.df[,1],invert=TRUE),]

cpgi.gr=GRanges(seqnames=cpgi.df[,1],
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))

#Creating a GRanges 
filePath=system.file("extdata",
                     "cpgi.hg19.chr21.bed",
                     package="compGenomRData")
filePath

df <- read.table(filePath)
df[grep("_", df[,1],invert = TRUE),]
df.gr = GRanges(seqnames = df[,1],
                ranges = IRanges(start=df[,2],
                                 end = df[,3]))


mcols(df.gr)=cbind(mcols(df.gr), DataFrame(name2=c("pax6","meis1","zic4")) )

#Useful methods to compare GRanges object
subsetByOverlaps(pk1.gr,cpgi.gr)
counts=countOverlaps(pk1.gr,cpgi.gr)
findOverlaps(pk1.gr,cpgi.gr)
n.ind=nearest(pk1.gr,cpgi.gr)
dists=distanceToNearest(pk1.gr,cpgi.gr,select="arbitrary")



#Creating a region to extract from bamfile corresponding to +/-1000 bp from start and end of promoter that are on chromosome 21
promoter.gr=tss.gr
start(promoter.gr)=start(promoter.gr)-1000
end(promoter.gr)  =end(promoter.gr)+1000
promoter.gr=promoter.gr[seqnames(promoter.gr)=="chr21"]



#Loading BAM file, setting region to scan and getting counts of region. Only loads region that was scanned. Requires a GRanges object to scan.
bamfilePath=system.file("extdata",
                        "wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam",
                        package="compGenomRData")

# get reads for regions of interest from the bam file
param <- ScanBamParam(which=promoter.gr)
counts=countBam(bamfilePath, param=param)


#Example GRanges object 
which.gr <- GRanges(c(
  + "seq1:1000-2000",
  + "seq2:100-1000",
  + "seq2:1000-2000"
  + ))
## equivalent:
## GRanges(
## seqnames = c("seq1", "seq2", "seq2"),
## ranges = IRanges(
## start = c(1000, 100, 1000),
## end = c(2000, 1000, 2000)
## )
#What to extract from the BAM file according to GRange object
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)


#Creating path to file and reading it into file according to parameters
bamfilePath=system.file("extdata", "ex1.bam", package="Rsamtools")
bam <- scanBam(bamFile, param=param)









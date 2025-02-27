library(GenomicRanges)
library(GenomeInfoDb)
library(GenomicAlignments)
library(ComplexHeatmap)
library(circlize)
library(rtracklayer)
library(Gviz)

data_path = system.file('extdata/chip-seq', package = 'compGenomRData')

chip_files = list.files(data_path, full.names=TRUE)


# fetch the chromosome lengths for the human genome
hg_chrs = getChromInfoFromUCSC('hg38')

# find the length of chromosome 21
hg_chrs = subset(hg_chrs, grepl('chr21$',chrom))

# downloaded hg_chrs is a data.frame object,
# we need to convert the data.frame into a named vector
seqlengths = with(hg_chrs, setNames(size, chrom))

# tileGenome function returns a list of GRanges of a given width, 
# spanning the whole chromosome
tilling_window = tileGenome(seqlengths, tilewidth=1000)

# unlist converts the list to one GRanges object
tilling_window = unlist(tilling_window)


# fetch bam files from the data folder
bam_files = list.files(
  path       = data_path, 
  full.names = TRUE, 
  pattern    = 'bam$'
)

# use summarizeOverlaps to count the reads
so = summarizeOverlaps(tilling_window, bam_files)

# extract the counts from the SummarizedExperiment
counts = assays(so)[[1]]

# calculate the cpm from the counts matrix
# the following command works because 
# R calculates everything by columns
cpm = t(t(counts)*(1000000/colSums(counts)))


# remove all tiles which do not contain reads
cpm = cpm[rowSums(cpm) > 0,]

# change the formatting of the column names
# remove the .chr21.bam suffix
colnames(cpm) = sub('.chr21.bam','',   colnames(cpm))

# remove the GM12878_hg38 prefix
colnames(cpm) = sub('GM12878_hg38_','',colnames(cpm))

# calculates the pearson correlation coefficient between the samples
correlation_matrix = cor(cpm, method='pearson')

#define the color palette which will be used in the heatmap
heatmap_col = circlize::colorRamp2(
  breaks = c(-1,0,1),
  colors = c('blue','white','red')
)

# plot the heatmap using the Heatmap function
Heatmap(
  matrix = correlation_matrix, 
  col    = heatmap_col
)


# list the bam files in the directory
# the '$' sign tells the pattern recognizer to omit bam.bai files
bam_files = list.files(
  path       = data_path, 
  full.names = TRUE, 
  pattern    = 'bam$'
)

# select the first bam file
chip_file = bam_files[1]


# load the genomic alignments package
library(GenomicAlignments)

# read the ChIP reads into R
reads = readGAlignments(chip_file)

# the reads need to be converted to a granges object
reads = granges(reads)


# extends the reads towards the 3' end
reads = resize(reads, width=200, fix='start')

# keeps only chromosome 21
reads = keepSeqlevels(reads, 'chr21', pruning.mode='coarse')


# convert the reads into a signal profile
cov = coverage(reads, width = seqlengths)

# change the file extension from .bam to .bigWig
output_file = sub('.bam','.bigWig', chip_file)

# export the bigWig output file
export.bw(cov, 'output_file')



# define the genome axis track
axis   = GenomeAxisTrack(
  range = GRanges('chr21', IRanges(1, width=seqlengths))
)

# convert the signal into genomic ranges and define the signal track
gcov   = as(cov, 'GRanges')
dtrack = DataTrack(gcov, name = "CTCF", type='l')

# define the track ordering
track_list = list(axis,dtrack)

# plot the list of browser tracks
# sizes argument defines the relative sizes of tracks
# background title defines the color for the track labels
plotTracks(
  trackList        = track_list, 
  sizes            = c(.1,1), 
  background.title = "black"
)


#Calculating plus-minus strand cross correlation using Jaccard index

# load the reads
reads = readGAlignments(chip_file)
reads = granges(reads)

# keep only the starting position of each read
reads = resize(reads, width=1, fix='start')

reads = keepSeqlevels(reads, 'chr21', pruning.mode='coarse')

calculate the coverage profile for plus and minus strand
reads = split(reads, strand(reads))

# coverage(x, width = seqlengths)[[1]] > 0 
# calculates the coverage and converts
# the coverage vector into a boolean
cov   = lapply(reads, function(x){
  coverage(x, width = seqlengths)[[1]] > 0
})

cov   = lapply(cov, as.vector)


# defines the shift range
wsize = 1:400

# defines the jaccard similarity
jaccard = function(x,y)sum((x & y)) / sum((x | y))

# shifts the + vector by 1 - 400 nucleotides and 
# calculates the correlation coefficient
cc = shiftApply(
  SHIFT = wsize, 
  X     = cov[['+']], 
  Y     = cov[['-']], 
  FUN   = jaccard
)

# converts the results into a data frame
cc = data.frame(fragment_size = wsize, cross_correlation = cc)


library(ggplot2)
ggplot(data = cc, aes(fragment_size, cross_correlation)) +
  geom_point() +
  geom_vline(xintercept = which.max(cc$cross_correlation), 
             size=2, color='red', linetype=2) +
  theme_bw() +
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('Shift in base pairs') +
  ylab('Jaccard similarity') 



#GC bias 

# fetches the chromosome lengths and constructs the tiles
library(GenomeInfoDb)
library(GenomicRanges)

hg_chrs        = getChromInfoFromUCSC('hg38')
hg_chrs        = subset(hg_chrs, grepl('chr21$',chrom))
seqlengths     = with(hg_chrs, setNames(size, chrom))

# tileGenome produces a list per chromosome
# unlist combines the elemenents of the list 
# into one GRanges object
tilling_window = unlist(tileGenome(
  seqlengths = seqlengths, 
  tilewidth  = 1000
))

# loads the human genome sequence
library(BSgenome.Hsapiens.UCSC.hg38)

# extracts the sequence from the human genome
seq = getSeq(BSgenome.Hsapiens.UCSC.hg38, tilling_window)


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


ggplot(
  data = gc, 
  aes(
    x = GC, 
    y = GM12878_hg38_CTCF_r1.chr21.bam
  )) +
  geom_point(size=2, alpha=.3) +
  theme_bw() +
  theme(
    axis.text  = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('GC content in one kilobase windows') +
  ylab('log10( cpm + 1 )') +
  ggtitle('CTCF Replicate 1')

# load the tidyr package
library(tidyr)

# pivot_longer converts a fat data.frame into a tall data.frame, 
# which is the format used by the ggplot package
gcd = pivot_longer(
  data      = gc, 
  cols      = -GC,
  names_to  = 'experiment',
  values_to = 'cpm'
)

# we select the ChIP files corresponding to the ctcf experiment
gcd = subset(gcd, grepl('CTCF', experiment))

# remove the chr21 suffix
gcd$experiment = sub('chr21.','',gcd$experiment)

ggplot(data = gcd, aes(GC, log10(cpm+1))) +
  geom_point(size=2, alpha=.05) +
  theme_bw() +
  facet_wrap(~experiment, nrow=1)+
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('GC content in one kilobase windows') +
  ylab('log10( cpm + 1 )') +
  ggtitle('CTCF Replicates 1 and 2')



# load the AnnotationHub package
library(AnnotationHub)

# connect to the hub object
hub = AnnotationHub()


# query the hub for the human annotation
AnnotationHub::query(
  x       = hub, 
  pattern = c('ENSEMBL','Homo','GRCh38','chr','gtf')
)

# retrieve the human gene annotation
gtf = hub[['AH61126']]

# extract ensemel chromosome names
ensembl_seqlevels = seqlevels(gtf)

# paste the chr prefix to the chromosome names
ucsc_seqlevels    = paste0('chr', ensembl_seqlevels)

# replace ensembl with ucsc chromosome names
seqlevels(gtf, pruning.mode='coarse') = ucsc_seqlevels

gtf = gtf[seqnames(gtf) == 'chr21']


# construct a GRangesList with human annotation
annotation_list = GRangesList(
  
  # promoters function extends the gtf around the TSS 
  # by an upstream and downstream amounts
  tss    = promoters(
    x = subset(gtf, type=='gene'), 
    upstream   = 1000, 
    downstream = 1000),
  exon   = subset(gtf, type=='exon'),
  intron = subset(gtf, type=='gene')
)


annotateReads = function(bam_file, annotation_list){
  
  library(dplyr)
  message(basename(bam_file))
  
  # load the reads into R
  bam    = readGAlignments(bam_file)
  
  # find overlaps between reads and annotation
  result = as.data.frame(
    findOverlaps(bam, annotation_list)
  ) 
  
  # appends to the annotation index the corresponding
  # annotation name
  annotation_name   = names(annotation_list)[result$subjectHits]
  result$annotation = annotation_name
  
  # order the overlaps based on the hierarchy
  result = result[order(result$subjectHits),]
  
  # select only one category per read
  result = subset(result, !duplicated(queryHits))
  
  # count the number of reads in each category
  # group the result data frame by the corresponding category
  result = group_by(.data=result, annotation)
  
  # count the number of reads in each category
  result = summarise(.data = result, counts = length(annotation))
  
  # classify all reads which are outside of 
  # the annotation as intergenic
  result = rbind(
    result, 
    data.frame(
      annotation = 'intergenic',
      counts     = length(bam) - sum(result$counts)
    )
  )
  
  # calculate the frequency
  result$frequency  = with(result, round(counts/sum(counts),2))
  
  # append the experiment name
  result$experiment = basename(bam_file)
  
  return(result)
}


# list all bam files in the folder
bam_files   = list.files(data_path, full.names=TRUE, pattern='bam$')

# calculate the read distribution for every file
annot_reads_list = lapply(bam_files, function(x){
  annotateReads(
    bam_file        = x, 
    annotation_list = annotation_list
  )
})

# collapse the per-file read distributions into one data.frame
annot_reads_df = dplyr::bind_rows(annot_reads_list)

# format the experiment names
experiment_name = annot_reads_df$experiment
experiment_name = sub('.chr21.bam','', experiment_name)
experiment_name = sub('GM12878_hg38_','',experiment_name)
annot_reads_df$experiment = experiment_name


ggplot(data = annot_reads_df, 
       aes(
         x    = experiment, 
         y    = frequency, 
         fill = annotation
       )) +
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_brewer(palette='Set2') +
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('Sample') +
  ylab('Percentage of reads') +
  ggtitle('Percentage of reads in annotation')
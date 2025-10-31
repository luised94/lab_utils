# Project development
## Intro
Just outlines what I am consider to complete in the project. TODO is things I am considering or currently working on. COMPLETE is tagged by date.

## TODO
### CODE focused
Single vs Paired end Analysis: Make scripts deal with the fact that paired end reads requires different submission requirements to different command line programs.

### Biology focused
Cross-correlation Analysis: Using phantompeakqualtools for strand cross-correlation (NSC/RSC metrics)
Peak Statistics: Using bedtools to analyze peak width distributions, distances between peaks, and genomic feature overlaps
Coverage Comparisons: Using deepTools multiBigwigSummary and plotCorrelation to compare signal profiles
Enrichment Analysis: Using GREAT or similar tools for genomic region enrichment
Motif Analysis: Using MEME-ChIP or HOMER for motif discovery in peaks
Peak Conservation: Analyzing conservation scores within peaks using phyloP/phastCons

## COMPLETE
### 2025-10-20
Was able to recapitulate most of Eaton 2010 algorithm for peak calling (280 peaks vs 267 peaks without taking the intersection of both replicates).

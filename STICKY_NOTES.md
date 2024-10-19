# STICKY NOTES
This note is just a section that documents thoughts and comments while I am developing the code in this directory.
I have separated it into dates.

## 2024 04 21
- Most of the bash code I wrote for my CHIP-seq analysis previously served only to run R code.
- Need to see how I read feature files.
- Can definitely the quality control with bash (deeptools or other)
- Think I finally understand how I was coding. I was essentially debugging and testing but leaving it in that state.
- Furthermore, I think I was using weird unintuitive approaches and building on top of it if it worked.
- Need to output environment info after I finish analysis with git info.
- Installed all packages in 1-package-installation.r. I am reworking this file so I will just have to remember to document installed and used pakcages.
- Was able to replace most the initial part of the pipeline with bash instead of R.
- Copied rsync files from previous code collecition  because it had useful parts to extract.
- Finally went through old code collection for bioinformatics analysis.

## 2024 04 22
- Was using assign for variable assignment in my previous code assignment, contributed to uninterpretability even though it is a way to assign variables programmatically
- Will implement a for loop or lapply probably for the track assignment and plotting, easy to understand
- Think I need to rerun the alignments because there may have been an ordering problem. Not sure since variable naming is done on basis of the file. 

## 2024 04 23
SPB comments
Overnight, Temp, Amount of antibody, same conditions as Kate from winston lab (zotero: @miller_winston23)
Dont think I have to repeat ORC two more times but probably will have to... 

## 2024 04 27
- Had some trouble with the text processing and used awk to make more robust. 
- Could come up with script that lets me know if the assumptions I make in my script are met. On the other hand, could just use same name instead of processing.
- Definitely thinking a little about design up front is worth it pero esto vale cuando haga otro proyecto. 

## 2024 04 28
- Creating a new branch for flow_cytomery since I am also attempting to use renv for package management.
Step-by-step: for future reference
'''{bash}
git branch feature/flow_cytometry

#Procede with feature development. In this case, creating the module for analysing flow_cytometry data.
git checkout feature/flow_cytometry
# Return to main
git checkout main 

'''

## 2024-07-11

Realized that I needed to edit the the OUT and ERR file creation in the next-generation-sequencign module in multiple files.

Used a nice sed command to check then in-line edit on the files.
Add the -i flag to perform inline editing.
Then I wanted to add certain lines after a particular line and at the end of the file.
Combined for loop with globbing and a multi-line sed command.

## 20240826
Have a decent set of features that I can use to compare to my data. Also have a good amount of reference to use for further comparisons and factor analysis.
Need to update the genome track plotting script and create all of the downstream analysis.
Still need to do quality control. 
Must work faster.

## 2024-09-06
Need to update the Gviz usage to ggplot2.
Considering using tidyverse and other packages for more readable and manageable code base.
Implement experiment definition in sampleConfig.R
Implement peak calling, motif analysis and comparisons and analysis.

## 2024-09-08
Output plots with tracks overlayed but with colors and different comparisons.
Update sampleConfig to ConfigAndExpr and main/function
Implement peak calling, motif analysis and comparisons and analysis.
Update track plotting to ggplot2.

## 2024-09-11
try lapply with all_tracks to highlight track and plot to see if any or all of the tracks trigger the error.
start new thread to get fresh help.

## 2024-09-12
Need to the unique labeling scheme for the comparison tracks. Task_Work_R_020: Update genomicTracks and implementation of unique name labeling
Think I can fix highlighting by making sure I use the right chromosome name.
Update the color probably.

## 2024-09-16
Updating the reference genome downloads file to use datasets.
The datasets command does not download the feature files for the W303 genome.
Wget triggers robots.txt.
Need to rework the indexing and the alignments but approach should be more robust.
See threads for potential solutions: 20240916_Task_Work_Bash_TBD: Update reference genome downloading using datasets cli. Revisit

## 2024-09-27
Have to update the generation of bigwig files to be only for S288C genome.
Need to see the plots for the generated bigwig files.
Need to adjust the pattern searching in 002_plotUserDefinedExperiments.R and other plotting functions. Will have to deal with this for other files such as fastqs and bam, etc.
Readjust the entire organization into functions.
See results using: find /home/luised94/data/240819Bel/logs -type f -name 202409270456*.out -exec vim {} +
Run 002_plotUserDefinedExperiments.R after bigwig files are generated.

## 2024-10-01
Create a comprehensive set of bigwig files by modifying the bamCompare and bamCoverage files to see the best way to visualize. 
Need to start working of peak calling.
Need to potentially add timeid or use distinct tags.

## 2024-10-18 
I think I should implement logging and then perform the reorganization into functions and scripts with updated logging.
Not super sure about the tests for each function or file but I guess that means I dont understand it enough.
Maybe quick runs using repl to get it to usable spot and then systematic test construction.

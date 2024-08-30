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

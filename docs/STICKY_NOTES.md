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
```{bash}
git branch feature/flow_cytometry

#Procede with feature development. In this case, creating the module for analysing flow_cytometry data.
git checkout feature/flow_cytometry
# Return to main
git checkout main
```

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

## 2024-10-23
Debating about where load and output type functions. I think load should go in file operations whole modify and output should go in table operations (or anything else that is being outputted.). This references the object that is being acted on. You load the file. You output the table.
If there is confusion about this because it requires, two things then it is likely that we should refactor it.
Will mark scripts files at the top to denote that I moved the functions to their respective files.
Hmm. It seems that separation of concerns also applies to testing.
Will further distinguish operations on filetype. For example, if file is tsv, fastq etc will be in tsv_operations, fastq_operations etc.
Need to break up functions in my config into functions and scripts and configuration.
May move xml creator into my config.
Didnt include the uninstallR script in the reorganization even tho I processed it to see. Not quite sure if I need it.
Reorganized till determineInputForArrayId.R

## 2024-10-24
Retain the FCS files but deal with them in the future.

## 2024-10-28
Currently have setup_bmc_experiment.R following a similar protocol to before. Initialize the samples in bmc_sample_grid_config.R. No longer output with comp_ columns. Instead they are found in the experiment config.
See 002__20241025_Task_Personal/Work_R/bash/Spaces/labutils_XXX: Using perplexity spaces to effectively update my configuration. Revisit. __project.
Search for current repository state one.
Have to add print statements for efficient debugging when of the bmc_sample_grid_config.

## 2024-10-29
Currently setup till transfer of fastq to bmc.
Mostly just run setup_bmc_experiment then transfer_bmc_experiment_directory_to_luria.
Havent granted authority to run as script. Run from repl but most source logging_utils.R and project_init.R, and project_config.R
Adjust bmc_sample_grid_config to setup the particular bmc experiment.
Then run  download_bmc_fastq_... .

## 2024-10-30
Remove and clean the directories 241007Bel and 241010Bel before redownloading.
Update run_alignment.
Clean up, then test slurm_wrapper with run_alignment.
Need to determine how to have log_file be empty to not mess up logging. Should be detected in the if statement at the bottom of log_message.
After testing download_bmc_fastq_to_user_bel_directory.sh, have to remove files referenced in current version: bmc_fastq_file_manager.sh,
Need to update BMC_CONFIG and its occurrences and calling of PROJECT_CONFIG.
I dont think verify_filesystem_path will work expected_fs is not equal to real_path.

## 2024-11-01
Consolidate the two core testing files.
Create the random lock files in bash and try to normalize the R and bash logging, lock and init.R.

## 2024-11-03
9641496
vim /home/luised94/logs/2024-11/bamcoverage/job_9641582/task_*/*.log
vim /home/luised94/logs/2024-11/bamcoverage/job_9641648/task_*/*.log

Need to doublecheck that samples are in the correct order than I submitted to the bmc.
Need to test that input sample name is being added correctly.

## 2024-11-07
Input samples could still be messed up. Ensure everything lines up during alignment, coverage calculations and plotting.

## 2024-11-11
Need to update the output directory to not nest plot directories. See script to plot all samples.

## 2024-12-15
Update the instances of:
- DEBUG_CONFIG enabled instances to single_file_mode.
- PLOT_CONFIG dimensions and name for particular specific scripts and rename t genome_track_config.

## 2024-12-16
Performed some benchmarking and simple scripts to analyze the 241010Bel dataset.
I tried normR and the following methods for isolating significant peaks:
1. taking the top 500.
2. taking the extreme peaks (qvals == 0)
3. using BH correction from the stats package.
4. trying to find the elbow of the q values.
None of these yielded peak sets that overlapped significantly with the eaton file.
The input file and the HM1108 looked noisy so changed to a cleaner input and a V5 sample and that also didnt work to produce biologically significant peak sets.
I plotted eaton and it looked ok. Furthremore, there is no input data so I cant use normR. Maybe I should have shifted at this time.
I think I need to plot particular regions, see if I can install macs2 manually, use another R package and finally complete the quality control.
I received the new repeat, so I can analyze that as well.
Need to reach out to an expert to see if they can help me.

## 2025-01-06
Came back from vacation to catastrophic tragedy. As alignment was occuring, I looked at logs and it seemed like the aignments were completing quickly.
After verification, a very small percentage of the reads were aligning to the reference genome.
After a discussion with shandon, she suggested I talk to the BMC to consult on possibilities.
I think I can still continue with the previous two repeats and check in with the BMC.
Need to ask about the numbered files as well. (245001-1 and 245001-2)
todo: Create the genomic tracks, finish the test module, update configs for the previous two repeats, prepare strategy and questions for BMC.

## 2025-04-03
I think I added a bunch of nonsense to most of my scripts with the loading modules and printing a lot of things, especially without having more control over what prints out to the screen.
They currently work at least but I may need to remove it.
Scripts to generate publication figures will be more targeted and streamlined either way.

## 2025-04-07
Similar problems with functions. They are probably doing too much but not as bad as last time.
Need to run the script for the bigwig files to see the effects of bam and bigwig processing. Should be fine to source from R after loading.

## 2025-04-23
Updated the pipeline completion scripts
Have a blacklist file now to start.
Need to rework the deduplication and shifting script to do blacklist removal after deduplication then do the shifting and fragment length prediction
Remove redundant operations from the deduplication script (like defining selected sample keys) by extracting metadata from the filename like I did in 003 and 004 scripts.
Remove or comment the logic for coverage calculation with blacklist

## 2025-04-25
Separated the file into different steps and tried to refactor most of the files.
This will definitely affect the downstream R scripts and will have to adjust how the metadata for each of the files is done. Could make the names independent instead of building them up.
E g name_of_bam_none_none.bam name_of_bam_deduped_none.bam etcetera.
Otherwise will have to handle missing columns. (maybe minimum is 3, "left align" the metadata then fill with default value, so like a lower part of a matrix from the diagonal match or fitted with the other half like a puzzle piece )
Need to find the plotting code but other than that I think its looking good.
Will comment the code to run as dry-run to double check then run.

## 2025-04-29
Almost done with the tests. Need to plot bigwig files and the peak file statistics
Bigwig file script runs, just need to plot and save to file to see if there are no errors.
For the peak files, can adapt the metadata processing but just need to go through all of the files to grab statistics and then decide on the plots.

## 2025-04-29
Finally got the data for the experiment.
I think I also got the pipeline completion scripts to run. Just need to see the output and generate the plots.
Considering creating one consolidated script that processes a file from fastq to bw, peaks, etcetera once I am sure I have the complete picture.
Need to find the plotting code for the test scripts.

## 2025-05-01
[DONE]TODO: Need to address the labels for the bigwig files. Either set manually in a mapping with the comparison quote (nested list with quote and the columns to use) or in the configuration file.
[IGNORE]TODO: Add border to the genome tracks?
[IGNORE]TODO: Create the map of the bad files?

## 2025-05-04
TODO: May need to do a big correlation analysis as a matrix heatmap. I will do that after looking at the visualizations and seeing if I need to rework some scripts.

## 2025-05-08
TODO: Need to adjust bamcoverage generation to see if file exists and exit gracefully if so.
```{troubleshooting}
Description: troubleshooting the plots for all of the dataframes being analyzed in analyze_peak_files.R
    > current_chrom_subset_df_processing_group_chromosome_heatmap_plot only has one column and has NA.
    > current_chrom_subset_df_processing_group_chromsome_counts_plot also doesnt show anything.
    > current_peak_subset_df_processing_group_2D_qvalue_vs_fold_enrich_plot Axis labels are weird. Move the processing group down, and the x-axis is too tight. Adjust the scale to be natural instead of 0,20. they dont reach 20. Otherwise, I think it is good.
    > Error, script exited.
        Plotting: current_peak_subset_df_processing_group_peak_num_vs_enrichment_plot
          Saving to: /home/luised94/data/preprocessing_test/plots/current_peak_subset_df_processing_group_peak_num_vs_enrichment_plot_Eaton.svg
        Error in stri_extract_first_regex(string, pattern, opts_regex = opts(pattern)) :
          Incorrectly nested parentheses in regex pattern. (U_REGEX_MISMATCHED_PAREN, context=`narrow|broad)`)In addition: There were 14 warnings (use warnings() to see them)
Conclusions:
    > Fixed the error with typical print debugging. Was subsetting using the peak logical vector instead of the chromosome.
    > Chromosome plots still very tight and hard to see if there are any differences. Have to remove the labels from the top and move them to the side.
    > Removing the 0-20 range helps. Need to adjust the labels.
    > Error had a missing parenthesis. Plot was clear tho.
```
```{troubleshooting}
Description: Troubleshooting peak and summary plots
    > Peak width distribution violin/box-plot looks ok but there really does not seem to be an effect.
Conclusions:
    > All plots work now. Could focus on cosmetics. Maybe I should output all of them and then do the cosmetics.
```

## 2025-05-14
QUESTION: Should I wrap the main logic around interactive? I have to add the cli vs interactive options.
QUESTION: Why does requireNamespace check take so long? Is it cached?

## 2025-05-14
TODO: Need to fix the bmc config and metadata to have the same column names

## 2025-05-20
Collect the plots side by side via illustrator then move on to flow cytometry.

## 2025-05-21
Something just looks wrong about the plots. They look terrible. I know duplication removal and normalization makes it look bad. Just want to see the raw data to make sure I am not missing anything. The terrible samples should remain the same.

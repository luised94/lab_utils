# Consolidating configuration values into a single experiment file
## Introduction
The repository has undergone several changes. I believe this is the last one that will involve a relatively dramatic change.
All _CONFIG variables will be moved to the experiment bmc config file.
This document contains the commands used to carry out this refactoring for future reference and learning.
Potential benefits:
- This should keep experiment configurations together, providing a single control point for each one.
- I will also not get a git notification everytime I update a script to change the settings, which I did not think about.
- Perform script operations based on experiment id provided.

## Operations
1. Centralize _CONFIG variables to template_bmc_config.R and delete the instances.
First, I move the _CONFIG files I found to the template_bmc_config.R file. I wanted to make sure I had moved all of them. I used quickfix list to create a list of locations that I could cycle through. After confirming, use a macro to remove the variables I wanted to move.
This was performed on files in the in the core_scripts/ directory only.
```{vim}
vimgrep /\w\+_CONFIG\s*<- list(/ `find . -type f -name "*.R" ! \( -name "template_bmc_config.R" -o -name "all_functions.R" -o -name "benchmark_peak_calling_normR.R" -o -name "peak_call_normR.R" -o -name "extract_bmcIDmetadata_process.R" -o -name "comparison_analysis.R" -o -name "functions_*" \)`
" Perform a deletion of a config object manually using <Ctrl+Q>f(%$d and record into macro.
cdo normal @q
" Cycle through rest of the entries to delete any missing entries (pattern doesnt necessarily capture all instances)
cn " Coupled with any manual deletions or additional @q runs.
```

The names used in the find command where added manually as I noticed what was in the quickfix list. Additional section delimeters where removed manually.
Confirm affected files using git status and git diff.
Still have leftover sections that I deleted manually using standard vim operations.
At this stage all of the scripts that I modified are broken and I will have to perform the operations on scripts in other branches, but good start.
Need to update the variable calls now, consolidate the configs, denest as much as possible.

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

2. Update, consolidate, rename, denest the config variables.
Carried out mainly with manual operations.
This will have downstream effects, requiring me to update all instances of the modified variable calls.
I wanted to have a set of statistics to reference before and after to see the changes quickly.
```{bash}
# Get all possible CONFIG references
grep -r "_CONFIG" . --include="*.R" > /tmp/all_configs.txt

# Extract unique variable names
cat /tmp/all_configs.txt | \
    perl -ne 'print "$1\n" if /([A-Z_]+_CONFIG)/' | \
    sort -u > /tmp/unique_configs.txt

# Show usage count for each
cat /tmp/all_configs.txt | \
    perl -ne 'print "$1\n" if /([A-Z_]+_CONFIG)/' | \
    sort | uniq -c | sort -nr
# Alternative with sed.
#grep -r "_CONFIG" . --include="*.R" | sed -n 's/.*\<\([A-Z][A-Z_]\+_CONFIG\)\>.*/\1/p' | sort -u
```

Now I would like to perform the renaming of the config variables in the files.
Open all files in nvim then create a quickfix list and a scratch buffer that will contain the lua script that will contain the script for renaming.

```{vim}
vimgrep /\w\+_CONFIG/j *.R
new
set buftype=nofile
" Use the following line to add the the unique variables that we are going to rename to save some time.
" read !grep -r "_CONFIG" . --include="*.R" | perl -ne 'print "$1\n" if /([A-Z_]+_CONFIG)/' | sort -u
```

```{lua}
-- Use Ctrl-Q 7j $A ="", to add the string at the end for 7 lines.
local config_renames = {
    DEBUG_CONFIG = "RUNTIME_CONFIG", -- More descriptive.
    EXPERIMENT_CONFIG = "EXPERIMENT_CONFIG",
    FASTQC_CONFIG = "FASTQC_CONFIG",
    GENOME_TRACK_CONFIG = "GENOME_TRACK_CONFIG",
    NORMR_CONFIG = "NORMR_CONFIG",
    PLOT_CONFIG = "GENOME_TRACK_CONFIG", -- More descriptive.
    TIME_CONFIG = "TIME_CONFIG",
    VIEWER_CONFIG = "VIEWER_CONFIG",
}
-- I encountered a few errors while trying to get the below to work. E5108, E486
-- I think the quickfix list wasnt initialized, buffers must be in memory, and if it is changed they must be written and saved.
local function process_configs()
    local qf_list = vim.fn.getqflist()
    if #qf_list == 0 then
        print("No quickfix entries found!")
        return
    end

    -- Save current buffer if modified
    if vim.bo.modified then
        vim.cmd('write')
    end

    for _, entry in ipairs(qf_list) do
        -- Force load buffer if not loaded, handle unsaved changes
        local filename = vim.fn.bufname(entry.bufnr)
        vim.cmd('silent! edit ' .. filename)

        if entry.bufnr and entry.lnum then
            local lines = vim.api.nvim_buf_get_lines(entry.bufnr, entry.lnum-1, entry.lnum, false)
            if #lines > 0 then
                local line = lines[1]
                for old_name, new_name in pairs(config_renames) do
                    if line:match(old_name) then
                        local cmd = string.format('buffer %d | %ds/%s/%s/gc',
                            entry.bufnr,
                            entry.lnum,
                            old_name,
                            new_name)
                        print("Executing: " .. cmd)
                        vim.cmd(cmd)
                        -- Save changes immediately
                        vim.cmd('update')
                        break
                    end
                end
            end
        end
    end
end

process_configs()
```

Rerun the grep commands without outputting to file to confirm that files were rename and see counts.
Was prepared to run git restore to rollback changes. Other options could have been to create an additional backup branch.
Now rename the variables in the template_bmc_config.R file.

3. Adding the template for argument parsing and validating presence of required configs.
Created the functions to make scripts process arguments.
Went to every script that needs to have the argument and config processing added and place #<SCRIPT_CONTROL> tag where I wanted to add the template block.
Create the quickfix list with either a find command or with the SCRIPT_CONTROL tag.
Create a new scratch buffer.
```{vim}
vimgrep /#<SCRIPT_CONTROL>/j ## " or a find command.
new
set buftype=nofile
```

Add the r code that we want to copy to the scratch buffer.

```{r}
# Source script control functions
source("functions_for_script_control.R")

# Parse arguments and validate configurations
args <- parse_args(commandArgs(trailingOnly = TRUE))
source(file.path("~/data", args[["experiment-id"]], "config/experiment_config.R"))
validate_configs(c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG"))
```

```{vim}
" Yank the entire thing into the plus register.
gg"+yG
" Then make a new buffer or in the same scratch buffer, record the macro.
new
" Go to SCRIPT_CONTROL, remove the line, paste the plus registers content, write.
qq/#<SCRIPT_CONTROL><CR>dd"+p:w<CR>q
" Use the quickfix list to perform the macro on all of the files.
cfdo normal @q
```

4. Update calls to CONFIG variables by creating a mapping in lua with quickfix list and template_bmc_config as reference.
```{vim}
" Open scratch buffer
new
set buftype=nofile
" Read in the variables that are referenced into the file.
r !grep -r "_CONFIG\\\$" . --include="*.R" | perl -ne 'print "$1\n" if /_CONFIG\$(\w+)/' | sort -u

" Insert the enclosing for the variable using lua syntax and add ="", to each line of the mapping.
" Create the quickfix list.
vimgrep /\w\+_CONFIG\$\w\+/j `find . -type f -name "*.R" ! \( -name "template_bmc_config.R" -o -name "all_functions.R" -o -name "benchmark_peak_calling_normR.R" -o -name "peak_call_normR.R" -o -name "extract_bmcIDmetadata_process.R" -o -name "comparison_analysis.R" -o -name "functions_*" \)`

" Go through each instance and set the mapping according to your new config variable name.
" Once you have set the mapping for the first time, you can remove all instances of mapping that you set.
g/RUNTIME_CONFIG\$debug_enabled/delete

" Complete the quickfix list.
" If there are any leftover you can remove those lines.
global/= "",/delete

" Recreate the quickfix list.
vimgrep /\w\+_CONFIG\$\w\+/j `find . -type f -name "*.R" ! \( -name "template_bmc_config.R" -o -name "all_functions.R" -o -name "benchmark_peak_calling_normR.R" -o -name "peak_call_normR.R" -o -name "extract_bmcIDmetadata_process.R" -o -name "comparison_analysis.R" -o -name "functions_*" \)`
```

Run the lua code again. It is interactive so each time you must confirm the change. Furthermore, it will try to correct all instances of the mapping in the line to your desired change. This will include, for example, function parameters such as verbose. You must skip this instances manually.

```{lua}

local config_renames = {
CATEGORIES = "CATEGORIES",
COLUMN_ORDER = "COLUMN_ORDER",
COMPARISONS = "COMPARISONS",
CONTROL_FACTORS = "CONTROL_FACTORS",
EXPERIMENTAL_CONDITIONS = "EXPERIMENTAL_CONDITIONS",
FILE_PATTERN = "FILE_PATTERN",
HEADER_PATTERN = "parse_header",
INVALID_COMBINATIONS = "INVALID_COMBINATIONS",
METADATA = "METADATA",
NORMALIZATION = "NORMALIZATION",
SAMPLE_CLASSIFICATIONS = "SAMPLE_CLASSIFICATIONS",
VERSION = "version_required",
VERSION_PATTERN = "version_pattern",
base_dir = "path_base",
chromosome = "process_chromosome",
comparison = "process_comparison",
display_time = "output_display_time",
dry_run = "output_dry_run",
enabled = "debug_enabled",
existing_version_limit = "version_max",
fastqc_pattern = "file_pattern",
files_to_process_idx = "process_file_index",
group = "process_group",
header_prefix = "parse_prefix",
height = "display_height",
interactive = "debug_interactive",
module_end = "parse_module_end",
module_names = "module_list",
module_reference_file = "path_module_ref",
module_separator = "parse_module_start",
output_suffix = "file_suffix",
patterns = "pattern_svg",
placeholder_suffix = "format_placeholder",
qc_subdir = "path_qc_dir",
samples_per_group = "process_samples_per_group",
save_plots = "output_save_plots",
single_file_mode = "process_single_file",
title_format = "title_dev_template",
track_name_format = "format_track",
validate_config = "debug_validate",
verbose = "debug_verbose",
width = "display_width",
}


local function process_configs()
    local qf_list = vim.fn.getqflist()
    if #qf_list == 0 then
        print("No quickfix entries found!")
        return
    end

    -- Save current buffer if modified
    if vim.bo.modified then
        vim.cmd('write')
    end

    for _, entry in ipairs(qf_list) do
        -- Force load buffer if not loaded, handle unsaved changes
        local filename = vim.fn.bufname(entry.bufnr)
        vim.cmd('silent! edit ' .. filename)

        if entry.bufnr and entry.lnum then
            local lines = vim.api.nvim_buf_get_lines(entry.bufnr, entry.lnum-1, entry.lnum, false)
            if #lines > 0 then
                local line = lines[1]
                for old_name, new_name in pairs(config_renames) do
                    if line:match(old_name) then
                        local cmd = string.format('buffer %d | %ds/%s/%s/gc',
                            entry.bufnr,
                            entry.lnum,
                            old_name,
                            new_name)
                        print("Executing: " .. cmd)
                        vim.cmd(cmd)
                        -- Save changes immediately
                        vim.cmd('update')
                        break
                    end
                end
            end
        end
    end
end

process_configs()
```

Now, before merging into main, I created a branch for main as backup before merging.
```{bash}
git checkout main
git checkout -b backup_main_250102_config_consolidation
```
This branch will contain the backup before merging the config_consolidation into main.

## Updating the config files for other experiments
I have updated the config file and done most of the updates I expected. Now I have to test the new module and that the configs are loaded and processed properly.

```{bash}
# Copied the documentation to the directories. These have the same samples and are repeats.
cp ~/data/241122Bel/documentation/* ~/data/241007Bel/documentation/

# Rename the files

# Directory path
DOC_DIR="$HOME/data/241007Bel/documentation"

# First, preview changes (dry run)
echo "Preview of changes to be made:"
find "$DOC_DIR" -type f -name "*241122Bel*" -exec bash -c '
    for file; do
        new_name="${file//241122Bel/241007Bel}"
        echo "Will rename:"
        echo "  FROM: $file"
        echo "  TO:   $new_name"
    done
' bash {} +

find "$DOC_DIR" -type f -name "*241122Bel*" -exec bash -c '
    for file; do
        new_name="${file//241122Bel/241007Bel}"
        mv -v "$file" "$new_name"
    done
' bash {} +
# Manually update the experiment id in the config files
nvim ~/data/*Bel/documentation/*bmc_config.R ~/lab_utils/core_scripts/template_bmc_config.R
```
Copy paste the config variables to each bmc_config.R.

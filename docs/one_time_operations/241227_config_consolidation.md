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

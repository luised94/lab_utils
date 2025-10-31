message("Sourcing .Rprofile...")
ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
#system("git rev-parse --show-toplevel", intern = TRUE)
#if command -v git &>/dev/null && git rev-parse --is-inside-work-tree &>/dev/null; then
#  git rev-parse --show-toplevel
#  return 0
#fi
if (!file.exists("~/.inputrc")){
  message("To set vi mode for the repl: ")
  message("Add set editing-mode vi to ~/.inputrc")
}
# Source renv when working on my local computer.
# Prevent renv activation while on cluster.
#nzchar(Sys.getenv("SLURM_JOBID"))
lscpu_output_chr <- system("lscpu", intern = TRUE)
cpu_info <- lscpu_output_chr[grep("Model name:", lscpu_output_chr)]
hostname <- system("hostname", intern = TRUE)
if(!grepl("Xeon", cpu_info) ||
   !nzchar(Sys.getenv("SLURM_JOBID")) &&
   hostname != "luria")
{
  source("renv/activate.R")
}

# Set to use nvim or vim as the viewer of the help page
options(help_type = "text")
options(pager = file.path(ROOT_DIRECTORY, "bash/helper/nvim_R_pager.sh"))

# Remove previously defined objects to avoid polluting environment.
rm(list = ls())

if(interactive()){
  edit_vim <- function(path){
    stopifnot(
      "File path does not exist." = file.exists(path)
    )
    vim_command <- paste("vim", path, sep = " ")
    system(vim_command)
  }
  git_pull <- function(){
    system("git pull")
  }
}

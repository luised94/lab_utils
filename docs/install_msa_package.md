---
Title: Installing R packages for sequence alignment and alignment visualization
Date: 2025-11-25
---

# Introduction
The purpose of this document is to write down the issues encountered while installing the msa package, the solutions attempted and if they succeeded.
The reason for installing the msa package is to generate alignments of fasta files and subsequently use these files to generate plots around regions of interest.

# First error: 2025-11-25
Reference: https://www.perplexity.ai/search/i-am-trying-to-install-msa-pac-sds7mLhxQbKfsITvkCom4w
## First attempt.
Try to add global environment setting for the flags.
```error
/usr/include/c++/11/bits/c++0x_warning.h:32:2: error: #error This file requires compiler and library support for the ISO C++ 2011 standard. This support must be enabled with the -std=c++11 or -std=gnu++11 compiler options.
   32 | #error This file requires compiler and library support
```

Solution:
```bash
mkdir -p ~/.R
nvim ~/.R/Makevars
# Write the line: CXX_STD = CXX11
```
```R
renv::install("msa") # Confirm proceed with Y
```

Result:
Did not work. I suspected it was because the environment was not loaded. I exited and restarted and checked some values.

## Second attempt
Set manually and install.
Solution:
```R
Sys.setenv(CXX_STD = "CXX11")
Sys.getenv("CXX_STD") # confirm
renv::install("msa") # Confirm proceed with Y
```

Result:
Did not work.

## Third attempt
Set via renv specifically. Setup a project specific R environ.
Solution:
```bash
nvim ./.Renviron
# Write CXX_STD=CXX11
```

Result:
Did not work.

## Fourth solution
Makevars with renv.
Solution:
Write to .Renviron.
    R_MAKEVARS_USER=renv/Makevars
Write to renv/Makevars.
    CXX_STD = CXX11
    CXX11 = g++ -std=c++11

Result:
Nope.

## Fifth solution
Modify the Makeconf being use by the compilation directly.
Find CXX11FLAGS and add the -std=c++11

## Sixth solution
Switch to DECIPHER package.

# Second error: 2025-11-25
Issues installing ggmsa. Depedency ggalt was removed in August...
Solution:
renv::install("https://cran.r-project.org/src/contrib/Archive/ggalt/ggalt_0.4.0.tar.gz")

That led to the error that I need libproj-dev libraries.
sudo apt-get install libproj-dev libgdal-dev

That led to the error:
E: dpkg was interrupted, you must manually run 'sudo dpkg --configure -a' to correct the problem.'

Ran the suggested command.
apt --fix-broken install

Then had to apt-get update.
Then installed the libproj-dev libraries.
Then retried renv for ggalt and ggmsa after that.

Very close to installing. Received warnng from treeio that could not resolve remote.
Received error for patchwork. object is_ggplot is not exported by ggplot2.

Install treeio from bioconductor.

# Remove potentially corrupt versions
remove.packages(c("patchwork", "ggplot2", "rlang"))

# Reinstall stable versions from CRAN
renv::install("rlang")
renv::install("ggplot2")
renv::install("patchwork")

Then retry: renv::install("ggmsa")

Got error for ggforce: Only available for C++14. Needed to adjust the Makeconf back.
Remove the modifications performed to install msa and rerun.
Finally, retry ggmsa.

Now got error with ggtree: object check_linewidth not found.
Reinstall ggplot2 with version 3.5.1.

Restart the session.

Still the same issue. Tried downgrading ggplot2.
Purge ggplot2 from cache.

Repurged ggplot2 and rlang, exited and retried.

Tried bypassing renv:
install.packages("BiocManager")
BiocManager::install("ggtree", lib = .libPaths()[1], force = TRUE)

Updated a bunch of packages.

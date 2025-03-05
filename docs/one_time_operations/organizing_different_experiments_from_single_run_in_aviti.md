---
title: Organizing multiple experiments from a single aviti run
---

Most of the experiments I perform are sequenced in different aviti runs for each experiment. In the case where multiple experiments are sequenced, they are returned in the same batch of sequencing results. This requires an additional step in processing.

```{bash}
# Typically, you just rsync according to the email from the bmc server.
# Substitute the variables with the respective values provided in email.
mkdir -p ~/data/${PROJECT_ID1}/fastq
mkdir -p ~/data/${PROJECT_ID2}/fastq
srun rsync -nav /net/${server}/data/bmc/public/${lab_name}/${PROJECT_ID} ~/data/${PROJECT_ID1}/fastq

# Based on how the samples are named, use globbing or find command to get the proper directories and then move them to the appropriate project id directory.
# Use echo to see dry-run"
echo mv ~/data/${PROJECT_ID1}/fastq/D25-1260{81..93}* ~/data/${PROJECT_ID2}/fastq/
```

Now the analysis can continue from the cleanup_bmc_directory.sh script.

# Experiments from same aviti run
- 250207Bel/250208Bel

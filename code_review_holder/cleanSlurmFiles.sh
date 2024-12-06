#!bin/bash

find . -maxdepth 1 -type f -name "slurm-*.out" -exec echo {} +

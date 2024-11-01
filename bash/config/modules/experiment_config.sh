#!/bin/bash

# File patterns and formats
export EXPERIMENT_FILE_PATTERN="*.sh"
export EXPERIMENT_NUMBER_FORMAT="%03d"
export DATE_FORMAT="%Y%m%d"

# File naming
export FILENAME_TEMPLATE="${DATE}_${EXPERIMENT_INDEX}_experiment.md"

# Paths
export TEMPLATE_DIR="templates"
export TEMPLATE_FILE="experiment_template.md"

#!/bin/bash
declare -A BMC_CONFIG=(
    [SOURCE_FS]="/net/bmc-pub17/data/bmc/public/Bell"
    [TARGET_FS]="/net/bmc-pub14/data"
    [MIN_SPACE_GB]=50
)
# BMC-specific settings
#[BMC_BASE_PATH]="/net/%s/data/bmc/public/Bell/%s"
#[BMC_FASTQ_DIR]="fastq"
#[BMC_DEFAULT_SERVER]="bmc-pub17"

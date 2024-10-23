#!/bin/bash

# Add to existing or create new configuration
declare -A DOC_CONFIG=(
    ["OUTPUT_FILE"]="repository_contents.xml"
    ["DEFAULT_VERBOSE"]=1
)

declare -A FILE_TYPES=(
    ["DOCUMENTATION"]="md"
    ["SCRIPTS"]="sh R lua py"
    ["CONFIG"]="xml json yaml yml"
)

declare -A EXCLUDED_PATTERNS=(
    ["VCS"]=".git"
    ["DEPRECATED"]="deprecatedCode"
    ["ENVIRONMENT"]="renv"
    ["TEMP"]="tmp temp"
)

declare -A XML_ENTITIES=(
    ["&"]="&amp;"
    ["<"]="&lt;"
    [">"]="&gt;"
    ["'"]="&apos;"
    ['"']="&quot;"
)

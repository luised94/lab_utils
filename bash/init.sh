#!/bin/bash
bash_function_directory="${HOME}/lab_utils/bash/functions/"

if [ ! -d $bash_function_directory ]; then
    echo "lab_utils function directory doesn't exist."
    echo "Use git to clone."
    exit 1
fi

while IFS= read -r -d $'\0' file; do
    echo "${file}"
    #source $file

done < <(find ${bash_function_directory} -type f -name "*.sh" -print0 )

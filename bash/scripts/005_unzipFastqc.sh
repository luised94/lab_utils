#!/bin/bash
#DESCRIPTION: Unzip the fastqc folders.
#USAGE: Doesnt work as a script but can be run from the command line.


#set -x

if [ $# -ne 1 ]; then
	echo -e "Provide directory to unzip fastqc files\n"
	exit 1
fi

DIRECTORY_TO_PROCESS="$(find -H $HOME/data -maxdepth 1 -type d -name "$1")"

if [ ! -d "$DIRECTORY_TO_PROCESS" ]; then
	echo -e "Error: Directory doesnt exist\n"
	exit 1
fi

echo "Will process $DIRECTORY_TO_PROCESS"
sleep 1 
find $DIRECTORY_TO_PROCESS/qualityControl -type f -name "*.zip" -exec echo {} \;

read -p "Proceed with unzipping these files? (y/n): " confirm && [[ $confirm == [yY] ]] || echo "Exiting" && exit 1

# If confirmed, proceed with deletion
echo "Unzipping files..."

#Can run this by itself on a particular directory. Not sure why it didnt work in the script.
find $DIRECTORY_TO_PROCESS/qualityControl -type f -name "*.zip" -exec sh -c 'unzip "{}"' \; 
#TODO Line above doesnt run inside script. Line below could be solution but havent tested. Could try moving directory, unzipping and moving back.
#find "$DIRECTORY_TO_PROCESS/qualityControl" -type f -name "*.zip" -exec sh -c '
#    cd "$(dirname "{}")" && unzip -o "$(basename "{}")" || echo "Failed to unzip: {}"' \;

echo "Files have been unzipped."

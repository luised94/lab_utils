mapfile -t FILENAMES < <(cat code.txt | rev | cut -d/ -f1 | rev | grep -v "#" | rev | sort | rev | uniq | grep -v "python3.10" | grep -v "renv" | grep -v ".Rproj")



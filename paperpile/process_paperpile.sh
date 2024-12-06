<<<<<<< HEAD
#STATUS: KEEP.
=======
<<<<<<< HEAD
#STATUS: KEEP.
=======
>>>>>>> 9eb69bf (feat(consolidate_into_main): Consolidate restructure into main to)
>>>>>>> quality_control
mapfile -t FILENAMES < <(cat code.txt | rev | cut -d/ -f1 | rev | grep -v "#" | rev | sort | rev | uniq | grep -v "python3.10" | grep -v "renv" | grep -v ".Rproj")



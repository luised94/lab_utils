#STATUS:
#!/bin/bash
#TODO need to create the script to rename files since the R_slurm or sh_node is not super useful. 
find . -maxdepth 1 -type f \( -name "*.sh" -o -name "*.R" \) | sed -E 's/_(sh|R)_(node|slurm)_/_/g'

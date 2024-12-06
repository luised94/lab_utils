#STATUS:
#Ran the commands manually as well instead of via script.
#Dry runs
find "$HOME/data/$1" -type d -name "*D24*" -exec echo "Deleting: {}" \; 
find "$HOME/data/$1" -type d -name "*D24*" | wc -l 


find "$HOME/data/$1" -type d -name "*D24*" -exec rm -rf {} + 
find "$HOME/data/$1" -type f -name "*unmapped*" -exec rm -rf + 


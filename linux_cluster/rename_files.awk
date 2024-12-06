#STATUS:
#BEGIN {
# Field separator designation.
    FS = "/"
}

{
# Filter the files that shouldnt be renamed.
# The other solution is to use a combination of find and grep commands to filter files and folders.
    if (index($NF, "_") == 0) next
    if ($NF ~ "rsync") next
    if ($0 ~ "test") {
        split($NF, parts, "_")
        new_name = parts[1] "_" "test" "_" parts[length(parts)]
    }
    else {
        split($NF, parts, "_")
        new_name = parts[1] "_" parts[length(parts)]
    }
    if (new_name != $NF) {
#    printf "Full name: %s\n", $0
#    printf "Old name: %s\n", $NF
# Use printf to create the command to pipe into sh.
# Have to escape all the quotes. Use substr to grab the path of the file and then combine with new name.
    printf "mv \"%s\" \"%s%s\"\n", $0, substr($0, 1, length($0)-length($NF)), new_name
    }
}

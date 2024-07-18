BEGIN {
    FS = "/"
}

{
    if (index($NF, "_") == 0) next
    if (index($NF, "rsync") == 1) next
    split($NF, parts, "_")
    new_name = parts[1] "_" parts[length(parts)]
    if (new_name != $NF) {
    printf "mv \"%s\" \"%s%s\"\n", $0, substr($0, 1, length($0)-length($NF)), new_name
    }
}

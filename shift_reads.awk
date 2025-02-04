#!/usr/bin/awk -f

BEGIN {
OFS = "\t"
#skipped = 0
    # Load chromosome sizes into AWK array
    while (getline < chrom_file){
            chrom_sizes[$1] = $2
    }
    close(chrom_file)
}

# Header lines
/^@/ { print; next }
# Skip reads with complex CIGAR operations
# $6 ~ /I|D/ { skipped++; next }

# Main processing
{
    chrom = $3
    orig_pos = $4
    strand = (and($2, 16) ? "reverse" : "forward")

    # Calculate new position
    new_pos = (strand == "forward") ? orig_pos + shift : orig_pos - shift

    # Clamping logic with debug checks
    if (chrom in chrom_sizes) {
        max_pos = chrom_sizes[chrom]
        if (new_pos < 1) new_pos = 1
        if (new_pos > max_pos) new_pos = max_pos
        clamped = (new_pos != orig_pos ? "CLAMPED" : "")
    } else {
        max_pos = "N/A"
        clamped = "NO_CHROM_DATA"
    }

    # Debug print every 100000th read
    if (NR % 100000 == 0) {
        printf "[AWK DEBUG] Read %d: %s:%d -> %s:%d (%s)\n",
            NR, chrom, orig_pos, chrom, new_pos, clamped > "/dev/stderr"
    }

    $4 = new_pos
    print
}

END {
        print "Skipped " skipped " reads with indels" > "/dev/stderr"
}

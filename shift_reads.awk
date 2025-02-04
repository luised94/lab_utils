#!/usr/bin/awk -f

BEGIN {
    OFS = "\t"
    # Initialize counters
    total_reads = 0
    forward_reads = 0
    reverse_reads = 0
    clamped_min = 0
    clamped_max = 0
    no_chrom_data = 0
    unmapped_reads = 0
    printf "SHIFT_STATS\t%s\n", strftime("%Y-%m-%d %H:%M:%S") >> log_file
    
    # Load chromosome sizes into AWK array
    while (getline < chrom_file){
        chrom_sizes[$1] = $2
    }
    close(chrom_file)
    
    # Initialize per-chromosome counters
    split("", reads_per_chrom)  # Empty associative array
}

# Header lines
/^@/ { print; next }

# Main processing
{
    total_reads++
    chrom = $3

    # Skip unmapped reads
    chrom == "*" { 
        unmapped_reads++
        # Output unmapped reads unchanged
        print
        next 
    }

    orig_pos = $4
    strand = (and($2, 16) ? "reverse" : "forward")
    # Track strand statistics
    if (strand == "forward") {
        forward_reads++
    } else {
        reverse_reads++
    }

    # Track per-chromosome counts
    reads_per_chrom[chrom]++

    # Calculate new position
    new_pos = (strand == "forward") ? orig_pos + shift : orig_pos - shift

    # Clamping logic with statistics
    if (chrom in chrom_sizes) {
        max_pos = chrom_sizes[chrom]
        if (new_pos < 1) {
            new_pos = 1
            clamped_min++
        }
        if (new_pos > max_pos) {
            new_pos = max_pos
            clamped_max++
        }
    } else {
        max_pos = "N/A"
        no_chrom_data++
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
    printf "TOTAL_READS\t%d\n", total_reads >> log_file
    printf "FORWARD_READS\t%d\t%.2f%%\n", forward_reads, (forward_reads/total_reads)*100 >> log_file
    printf "REVERSE_READS\t%d\t%.2f%%\n", reverse_reads, (reverse_reads/total_reads)*100 >> log_file
    printf "CLAMPED_MIN\t%d\t%.2f%%\n", clamped_min, (clamped_min/total_reads)*100 >> log_file
    printf "CLAMPED_MAX\t%d\t%.2f%%\n", clamped_max, (clamped_max/total_reads)*100 >> log_file
    printf "NO_CHROM_DATA\t%d\t%.2f%%\n", no_chrom_data, (no_chrom_data/total_reads)*100 >> log_file
    printf "UNMAPPED_READS\t%d\t%.2f%%\n", unmapped_reads, (unmapped_reads/total_reads)*100 >> log_file
    printf "PER_CHROMOSOME_STATS\n" >> log_file
    for (chrom in reads_per_chrom) {
        printf "CHROM\t%s\t%d\t%.2f%%\n", 
            chrom, reads_per_chrom[chrom], 
            (reads_per_chrom[chrom]/total_reads)*100 >> log_file
    }
    printf "END_STATS\n" >> log_file

}

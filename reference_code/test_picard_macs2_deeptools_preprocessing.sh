MACS2_ENV="macs2_env"
OUTDIR="macs2_test_results"
GENOME_SIZE="1.2e7"  # S. cerevisiae genome size
# From 241010Bel, any of the scripts to find appropriate files and set them manually. Should be 18 and 3. Input files were noisy and spiky relative to expected flat input.
TEST_SAMPLE="/home/luised94/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam"
INPUT_CONTROL="/home/luised94/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam"
# From 100303Bel, file is from 2010 Eaton reference paper.
# Samples are mismatched with the metadata order.
# Correct index for HM1108 antibody is 4.
REF_SAMPLE="/home/luised94/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
OUTPUT_PREFIX="macs2_test"

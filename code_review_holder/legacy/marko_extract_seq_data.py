import numpy as np
import pandas as pd
import pickle
import pysam

#define two helper functions to identify the base value (nucleotide) and the corresponding confidence (how sure are we that the value is right):
def most_common_base(lst):
    """
    Identifies the most common base in a given list of bases. Input is a list of bases. Output is most common base or 'N'
    """
    try:
        return max(set(lst), key=lst.count)
    except:
        return 'N'
def confidence(lst):
    """
    Calculates the confidence for most common base in a given list. conf = (most_common_base_count / total_read_count_for_base)*100.
    Input is list of bases. Output is confidence value for most common base or NaN.
    """
    try:
        return round(float(lst.count(most_common_base(lst)))/float(len(lst)), 2)*100
    except:
        return np.nan
#this function reads the sequencing files and extracts relevant informantion
def extract_seq_data(df, filename, suffix):
    '''
    Function extracts sequencing data from .bam or .sam files. Inputs are empty Pandas DataFrame with just column names; filename for the .bam/.sam file; and suffix
    for the biological type of data ('Case' vs 'Control' | 'Hit' vs 'Control' | or something like that, users joice)
    '''
    #list of Yeast chromosomes
    chromosomes = 'chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI'.split()
    #temporary lists to hold data
    chrom_name = []
    positions = []
    coverage = []
    base_confidence = []
    bases = []
    #file that holds sequencing data
    bamfile = pysam.AlignmentFile(filename, "rb" )
    for chromosome in chromosomes:
        print(chromosome, filename)
        #loop through the sequencing file and extracting relevant informantion:
        #pileupcolumn represents all the reads in the SAM/BAM file that map to a single base in the reference sequence.
        for pileupcolumn in bamfile.pileup(chromosome, truncate=True):
            chrom_name.append(chromosome) #chromosome name.
            positions.append(pileupcolumn.pos) #nucleotide position
            coverage.append(pileupcolumn.n) #coverage for the given base. how many times it was sequecned
            bases_in_reads = [] #hold all the seq info for a given base.
            #loop through all the reads that correspond to a given position
            for pileupread in pileupcolumn.pileups:
                #this loops through all the individual reads for a given base and makes a list of them.
                try:
                    bases_in_reads.append(pileupread.alignment.query_sequence[pileupread.query_position])
                except:
                    bases_in_reads.append('N')
            #uses the bases_in_reads list to find the likeliest nucleotide match for the base and the confidence value for that
            bases.append(most_common_base(bases_in_reads))
            base_confidence.append(confidence(bases_in_reads))

    df = pd.DataFrame(
        {'chromosome': chrom_name,
         'base_nr': positions,
         'base_value_{}'.format(suffix):bases,
         'coverage_{}'.format(suffix): coverage,
         'confidence_{}'.format(suffix):base_confidence
        })
    bamfile.close()
    return df

#pickle the DataFrames for all three experiments, as this part takes a long time...
for index in range(1,10):
    print('Analyzing experiment nr: {}'.format(index))
    #create an empty DataFrame to store _hit data and then call the function to fill the DataFrame with the data
    df_hit = pd.DataFrame({'chromosome': [], 'base_nr': [], 'base_value_hit': [], 'coverage_hit': [], 'confidence_hit': []})
    df_hit = extract_seq_data(df_hit, "/mnt/c/Users/Luis/Desktop/ORC4RA_Supp/Exp_{}_hit_sorted.bam".format(index), 'hit')
    #create an empty DataFrame to store _ctrl data and then call the function to fill the DataFrame with the data
    df_ctrl = pd.DataFrame({'chromosome': [], 'base_nr': [], 'base_value_ctrl': [], 'coverage_ctrl': [], 'confidence_ctrl': []})
    df_ctrl = extract_seq_data(df_hit, "/mnt/c/Users/Luis/Desktop/ORC4RA_Supp/Exp_{}_ctrl_sorted.bam".format(index), 'ctrl')
    #merge the two DataFrames, so we could better compare the data
    df_final = pd.merge(df_hit, df_ctrl)
    #let's pickle these so we would not need to run it 8 times:
    df_final.to_pickle('/mnt/c/Users/Luis/Desktop/ORC4RA_Supp/pickles/df_final_exp{}.pickle'.format(index))

#Exp_{}_hit.bam

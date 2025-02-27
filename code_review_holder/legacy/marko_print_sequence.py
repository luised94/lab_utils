import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

def print_sequence(experiment_nr, chromosome, start, end):
    print('experiment_nr:', experiment_nr)
    df_final = pd.read_pickle('/mnt/c/Users/Luis/Desktop/ORC4RA_Supp/pickles/df_final_exp{}.pickle'.format(experiment_nr))

    print(df_final.head())

    df_final['print_area'] = ((df_final['chromosome'] == chromosome) &
                            (df_final['base_nr'] > start) &
                            (df_final['base_nr'] < end))
    df_gene =  df_final[df_final['print_area']].copy()

    ctrl = list(df_gene['base_value_ctrl']) #gets sequence for area
    hit = list(df_gene['base_value_hit']) #gets sequence for area
    ctrl_sequence = ''
    hit_sequence = ''
    for i in range(len(ctrl)):
        ctrl_sequence += ctrl[i]
        hit_sequence += hit[i]

    return ctrl_sequence, hit_sequence

def print_coverage(experiment_nr, chromosome, start, end):
    print('experiment_nr:', experiment_nr)
    df_final = pd.read_pickle('/mnt/c/Users/Luis/Desktop/ORC4RA_Supp/pickles/df_final_exp{}.pickle'.format(experiment_nr))

    print(df_final.head())

    df_final['print_area'] = ((df_final['chromosome'] == chromosome) &
                            (df_final['base_nr'] > start) &
                            (df_final['base_nr'] < end))
    df_gene =  df_final[df_final['print_area']].copy()

    ctrl = list(df_gene['coverage_ctrl']) #gets coverage for area
    hit = list(df_gene['coverage_hit']) #gets coverage for area

    return ctrl, hit

action = str(input("enter 's' to print sequence; 'c' to print coverage: "))
experiment = int(input('please enter experiment nr from 1-9: '))
chromosome = str(input('please enter chromosome in format chrXVI: '))
start = int(input('please enter start nucleotide position as integer: '))
end = int(input('please enter end nucleotide position as integer: '))

if action == 's':
    ctrl_sequence, hit_sequence = print_sequence(experiment, chromosome, start, end)
    print('hit sequence:')
    print(hit_sequence + '\n')
    print('ctrl sequence:')
    print(ctrl_sequence)
elif action == 'c':
    ctrl_coverage, hit_coverage = print_coverage(experiment, chromosome, start, end)
    to_plot = input("Do you want to plot the coverage ('y' or 'n'): ")
    if to_plot == 'y':
        x_val = []
        count = start
        for a in ctrl_coverage:
            x_val.append(count)
            count += 1

        sns.scatterplot(x_val, ctrl_coverage)
        sns.scatterplot(x_val, hit_coverage)
        plt.legend(labels=['coverage_ctrl', 'coverage_hit'])
        plt.title("Readcount of experiment {}: ({}: {}-{})".format(experiment,chromosome, start, end))
        plt.xlabel("Chromosome: {} (bp)".format(chromosome))
        plt.ylabel("Readcount")
        plt.show()
    else:
        print(ctrl_coverage)
        print(hit_coverage)

else:
    print('please press either s or c')

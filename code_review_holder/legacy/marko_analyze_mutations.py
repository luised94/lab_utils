import pandas as pd
import pickle
#Modify variable to adjust for confidence filter
confidence_value = 75

def find_candidate_genes(experiment_nr):
    print('experiment_nr:', experiment_nr)
    df_final = pd.read_pickle('/mnt/c/Users/Luis/Desktop/ORC4RA_Supp/pickles/df_final_exp{}.pickle'.format(experiment_nr))

    df_final['is_different'] = ((df_final['coverage_hit'] > 10) & #we want to compare well sequenced areas of the genome (readcount > 10)
                            (df_final['confidence_hit'] > confidence_value) & #high confidence means that mutaion was present in most reads
                            (df_final['coverage_ctrl'] > 10) & #we want to compare well sequenced areas of the genome (readcount > 10)
                            (df_final['confidence_ctrl'] > confidence_value) & #high confidence means that mutaion was present in most reads
                            (df_final['base_value_hit'] != df_final['base_value_ctrl'])) #we expect the mutation to be different between "hit" and "ctrl"
    df_candidate_genes =  df_final[df_final['is_different']].copy()


    for index in df_candidate_genes.index:
        # df_candidate_genes.at[index, 'gene_name'] = row["sgdAlias"]
        # df_candidate_genes.at[index, 'gene_sequence'] = row["chromosome.sequence.residues"]
        df_candidate_genes.at[index, 'exp_nr'] = 'Exp_' + str(experiment_nr)
    print(df_candidate_genes.head())
    return df_candidate_genes

df_exp1_genes = find_candidate_genes(1)
df_exp2_genes = find_candidate_genes(2)
df_exp3_genes = find_candidate_genes(3)

df_exp4_genes = find_candidate_genes(4)
df_exp5_genes = find_candidate_genes(5)
df_exp6_genes = find_candidate_genes(6)

df_exp7_genes = find_candidate_genes(7)
df_exp8_genes = find_candidate_genes(8)
df_exp9_genes = find_candidate_genes(9)
#
# df_all_candidate_genes = pd.concat([df_exp1_genes, df_exp2_genes, ])

df_all_candidate_genes = pd.concat([df_exp1_genes, df_exp2_genes, df_exp3_genes,
                                    df_exp4_genes, df_exp5_genes, df_exp6_genes,
                                    df_exp7_genes, df_exp8_genes, df_exp9_genes, ])
df_all_candidate_genes.to_csv('all_candidate_genes_'+str(confidence_value)+'.csv')

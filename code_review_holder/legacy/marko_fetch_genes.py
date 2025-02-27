from __future__ import print_function
from intermine.webservice import Service
import pandas as pd
import numpy as np

#Modify variable to adjust for file to use
confidence_value = 75
df = pd.read_csv('/mnt/c/Users/Luis/Desktop/ORC4RA_Supp/pickles/all_candidate_genes_'+str(confidence_value)+'.csv', index_col='index')
df["gene_name"] = 'LOOK UP'
df["gene_name2"] = 'LOOK UP'
df["gene_name3"] = 'LOOK UP'
print(df.head())

for index, row in df.iterrows():
        service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
        query = service.new_query("Gene")
        query.add_view(
            "chromosome.primaryIdentifier", "primaryIdentifier", "secondaryIdentifier",
            "featureType", "symbol", "name", "sgdAlias", "organism.shortName",
            "qualifier", "chromosomeLocation.start", "chromosomeLocation.end",
            "chromosomeLocation.strand", "description"
        )
        query.add_sort_order("Gene.primaryIdentifier", "ASC")
        query.add_constraint("chromosome.primaryIdentifier", "=", row['chromosome'])
        query.add_constraint("chromosomeLocation.start", "<=", str(row['base_nr']))
        query.add_constraint("chromosomeLocation.end", ">", str(row['base_nr']))

        for gene in query.rows():
            # row['gene_name'] = gene['sgdAlias']
            df.loc[index,'gene_name'] = gene['sgdAlias']
            df.loc[index,'gene_name2'] = gene['name']
            df.loc[index,'gene_name3'] = gene['symbol']
            print(gene['sgdAlias'])

df.to_csv('genes_fetch_'+str(confidence_value)+'.csv')
print(df.tail())

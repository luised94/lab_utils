#This script outputs 
#Unification of Protein Abundance Datasets Yields a Quantitative Saccharomyces cerevisiae Proteome.
import os

os.chdir("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/")

from query_functions import query_to_df
from datetime import date
# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "proteins.symbol", "secondaryIdentifier", "chromosome.primaryIdentifier",
    "proteins.median",  "proteins.MAD", "proteins.units"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "B")
query.add_constraint("Gene", "IN", "ALL_Verified_Uncharacterized_Dubious_ORFs", code = "A")


filename =str(date.today()) + '_Protein_Median_Abundance.xlsx'
column_names = ['Protein_Name','Systematic_Name','Chromosome','Median_Abundance','Absolute Median Deviation','Units']
views_list = ["proteins.symbol", "secondaryIdentifier", "chromosome.primaryIdentifier",
    "proteins.median",  "proteins.MAD", "proteins.units"] 


#query_to_df(query, views_list, column_names, to_excel=1, Filename=str(date.today())+'query.xlsx'):
df = query_to_df(query, views_list, column_names, Filename=filename, To_excel=0) 
file=open(filename, 'w')
df.to_excel(filename, index=False)
file.close()
# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("B and A")

# for row in query.rows():
#     print(row["primaryIdentifier"], row["proteins.symbol"], row["secondaryIdentifier"], \
#         row["proteins.median"], row["proteins.units"], row["proteins.MAD"], \
#         row["proteins.proteinAbundance.publication.pubMedId"])


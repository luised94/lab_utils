
"""
Created on Sun Apr 26 12:29:57 2020

@author: liusm

Examples of Queries from yeastmine

"""


# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

"""
Use to import functions for queries

"""
# import os

# file = os.chdir(r"C:\Users\liusm\Desktop\ARS_and_Gene_Analysis\YeastMine")

# from query_functions import query_to_df

"""

Yeast Intermine Service

"""
# The following two lines will be needed in every python script:
from datetime import date
from intermine.webservice import Service

yeast_url = "https://yeastmine.yeastgenome.org/yeastmine/service"
service = Service(yeast_url)
service = Service(yeast_url, username = "luised94@mit.edu", password = "Jkiu56rt!")

"""
Example for Gene Query

"""
# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "secondaryIdentifier", "symbol", "name", "sgdAlias", "phenotypeSummary",
    "chromosome.primaryIdentifier", "chromosomeLocation.start",
    "chromosomeLocation.end", "chromosomeLocation.strand"
)
#"upstreamIntergenicRegion.chromosomeLocation.start","upstreamIntergenicRegion.chromosomeLocation.end",upstreamIntergenicRegion.length"
    
# This query's custom sort order is specified below:
query.add_sort_order("chromosome.primaryIdentifier")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "A")
query.add_constraint("Gene", "IN", "Uncharacterized_Verified_ORFs", code = "B") 
query.add_constraint("primaryIdentifier", "!=", "chrmt", code = "C")#can change "Verified_ORFs" parameter to include dubious ORFs
# Uncomment and edit the code below to specify your own custom logic: Default is AND between constraints
# query.set_logic("A")

"""
Example for ARS Query

"""

query = service.new_query("ARS")


# The view specifies the output columns. Add view using syntax from yeast intermine and place strings in order
#you want column output to be. 
#ARS in yeastmine dont have directionality. Dont add "chromosomeLocation.strand"
query.add_view(
     "primaryIdentifier", "secondaryIdentifier", "symbol",
    "chromosome.primaryIdentifier", "chromosomeLocation.start",
    "chromosomeLocation.end", "sequence.residues"
)

# This query's custom sort order is specified below:
query.add_sort_order("chromosome.primaryIdentifier")

# You can edit the constraint values below
query.add_constraint("ARS", "IN", "ARSs", code = "A")

"""
Example for Chromosome Query

"""

query = service.new_query("Chromosome")


query.add_view("primaryIdentifier", "length")

query.add_sort_order("primaryIdentifier")

query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "A")
query.add_constraint("primaryIdentifier", "!=", "chrmt", code = "B")


"""
Example for DubORFs Query

"""

query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "secondaryIdentifier", "symbol", "name", "sgdAlias", "phenotypeSummary",
    "chromosome.primaryIdentifier", "chromosomeLocation.start",
    "chromosomeLocation.end", "chromosomeLocation.strand"
)
#"upstreamIntergenicRegion.chromosomeLocation.start","upstreamIntergenicRegion.chromosomeLocation.end",upstreamIntergenicRegion.length"
    
# This query's custom sort order is specified below:
query.add_sort_order("chromosome.primaryIdentifier")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "A")
query.add_constraint("Gene", "IN", "Dubious_ORFs", code = "B") 
query.add_constraint("primaryIdentifier", "!=", "chrmt", code = "C")#can change "Verified_ORFs" parameter to include dubious ORFs



query = service.new_query("Gene")
query.add_view(
    "secondaryIdentifier", "symbol", "proteins.sequence.residues",
    "proteins.sequence.length", "proteins.symbol", "locations.start",
    "locations.end", "chromosome.primaryIdentifier"
)
query.add_constraint("featureType", "=", "ORF", code = "A")
query.add_constraint("featureType", "=", "pseudogene", code = "C")
query.add_constraint("featureType", "=", "transposable_element_gene", code = "D")
query.add_constraint("Gene", "IN", "Verified_ORFs", code = "B")
query.set_logic("B and (A or C or D)")


"""
Outputting Query as CSV 
Adjust variables to query of interest

"""
query_name = "Chr_info_"

today = str(date.today())
filename =  query_name + today +'.csv'
f = open(filename, 'w')
views_list = ["primaryIdentifier", "length"]
header_string = 'Chr,Length\n'
f.write(header_string)

for row in query.rows():
    count = 1
    row_string = ''
    for view in views_list:
        string_holder = str(row[view]).replace(',', '_' ) #string holder to be able to modify twice, replace , with _ to write to csv
        if count == len(views_list):
            row_string += string_holder.replace('\n', '') + '\n'  #replace all new line characters in output of rows, particularly phenotype summary 
        else:
            row_string += string_holder.replace('\n', '') + ','
        count += 1
    f.write(row_string)
        
f.close()

""" 

Outputting Query as DF or CSV through DF (Requires importing function query_to_df )
Adjust variables to proper query you want to output

"""
column_names = ["Identified", "Lenghth"]
to_csv = 0 
query_to_df(query, views_list, column_names, to_csv)
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 14:04:00 2020

@author: liusm
"""

import os

os.chdir("C:/Users/liusm/Desktop/ARS_and_Gene_Analysis/YeastMine/")

from query_functions import query_to_df
# os.chdir("Path to dir that you want if you want to change back")

from datetime import date
from intermine.webservice import Service

yeast_url = "https://yeastmine.yeastgenome.org/yeastmine/service"
service = Service(yeast_url)
service = Service(yeast_url, username = "luised94@mit.edu", password = "Jkiu56rt!")

query = service.new_query("ARS")


# The view specifies the output columns. Add view using syntax from yeast intermine and place strings in order
#you want column output to be. 
#ARS in yeastmine dont have directionality. Dont add "chromosomeLocation.strand"
query = service.new_query("ARS")

# The view specifies the output columns
query.add_view(
    "primaryIdentifier", "secondaryIdentifier", 
    "chromosome.primaryIdentifier", "chromosomeLocation.start",
    "chromosomeLocation.end", "sequence.length", "sequence.residues"
)

# "downstreamIntergenicRegion.chromosomeLocation.start",
#     "downstreamIntergenicRegion.chromosomeLocation.end",
#     "upstreamIntergenicRegion.chromosomeLocation.start",
#     "upstreamIntergenicRegion.chromosomeLocation.end"
# You can edit the constraint values below
query.add_constraint("ARS", "IN", "ARSs", code = "A")

# This query's custom sort order is specified below:
query.add_sort_order("chromosome.primaryIdentifier", "ASC")
query.add_sort_order("chromosomeLocation.start", "ASC")



#Info for file details and output
today = str(date.today())
filename = "ARS"+ "_" +'fasta_' + today +'.fasta'
f = open(filename, 'w')

view_for_fasta_heading = ["secondaryIdentifier", "chromosome.primaryIdentifier", "chromosomeLocation.start", "chromosomeLocation.end"]
view_for_fasta_content  = "sequence.residues"
short_Description = "ARS Sequence"

#Writing output to fasta
for row in query.rows():
    row_string = '>'
    for view in view_for_fasta_heading:
        row_string += str(row[view]) + " "        
    f.write(row_string + short_Description)
    
    sequence_to_write = row[view_for_fasta_content]
    
    index_to_add_new_line = [x for x in range(len(sequence_to_write)+1) if x % 40 == 0]
    
    for index in index_to_add_new_line:
        sequence_to_write = sequence_to_write[:index] + '\n' + sequence_to_write[index:]
        
    f.write(sequence_to_write + "\n")

f.close()
#query_to_df(query, views_list, column_names, to_csv):
views_list = ["primaryIdentifier", "secondaryIdentifier", 
    "chromosome.primaryIdentifier", "chromosomeLocation.start",
    "chromosomeLocation.end", "sequence.length", "sequence.residues"]
column_names = ["Sys_Name", "Standard_Name", "Chr", "Start", "End", "Length", "Sequence"]

ARS_df = query_to_df(query, views_list, column_names, 0)
query_to_df(query, views_list, column_names, 1)

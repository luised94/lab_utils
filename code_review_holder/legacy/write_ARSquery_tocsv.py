# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 15:09:22 2020

@author: liusm
"""
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/


# The following lines will be needed in every python script:
from datetime import date
from intermine.webservice import Service

yeast_url = "https://yeastmine.yeastgenome.org/yeastmine/service"
service = Service(yeast_url)


# Get a new query on the class (table) you will be querying:
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

#Create list to iterate through when writing to csv, update when you update views.
views_list = [ "primaryIdentifier", "secondaryIdentifier", "symbol",
    "chromosome.primaryIdentifier", "chromosomeLocation.start",
    "chromosomeLocation.end",  "sequence.residues"]


#Create and open file to write query to csv
today = str(date.today())
filename = 'ARS_info_' + today +'.csv'
f = open(filename, 'w')

#Create header string with column names. Has to be changed depending on columns and prefered names for them
header_string = 'DBID,SystematicName,StandardName,ChrNum,Genome_Start,Genome_end,Sequence\n'
f.write(header_string)
    
#Writing out the row data as csv. For every row, add the data as strings with comma (,) at end unless
#it is last element, in which case add new line to be able to add other row data below. Write to file. 
for row in query.rows():
    count = 1
    row_string = ''
    for view in views_list:
        if count == len(views_list):
            row_string += str(row[view]).replace(',', '/') + '\n'
        else:
            row_string += str(row[view]).replace(',', '/') + ','
        count += 1
    f.write(row_string)
        
f.close()

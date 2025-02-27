# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 20:10:49 2020

@author: liusm
"""

# from datetime import date
from intermine.webservice import Service

from datetime import date
from query_functions import query_to_df

service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
query = service.new_query("Chromosome")

views = ["primaryIdentifier", "length"]
column_names = ['Chr','Length']
query.add_view("primaryIdentifier", "length")

query.add_sort_order("primaryIdentifier")

query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "A")
query.add_constraint("primaryIdentifier", "!=", "chrmt", code = "B")


    
    
    
chrom_dataframe = query_to_df(query, views,column_names,0)
 
today = str(date.today())
filename = 'Chr_info_' + today +'.csv'
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

    
    
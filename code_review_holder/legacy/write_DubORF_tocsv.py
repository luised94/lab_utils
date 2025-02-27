# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:27:00 2020

@author: liusm
"""

# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/


# The following two lines will be needed in every python script:
from datetime import date
from intermine.webservice import Service

yeast_url = "https://yeastmine.yeastgenome.org/yeastmine/service"
service = Service(yeast_url)

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
query.add_constraint("Gene", "IN", "Dubious_ORFs", code = "B") 
query.add_constraint("primaryIdentifier", "!=", "chrmt", code = "C")#can change "Verified_ORFs" parameter to include dubious ORFs
# Uncomment and edit the code below to specify your own custom logic: Default is AND between constraints
# query.set_logic("A")

# for row in query.rows():
#     print(row["secondaryIdentifier"], row["symbol"], row["name"], row["sgdAlias"], \
#         row["chromosome.primaryIdentifier"], \
#         row["chromosomeLocation.start"], row["chromosomeLocation.end"], \
#         row["chromosomeLocation.strand"])
 
#Make list to be able to iterate through row data of query. To modify, add view in order you want    
views_list = ["secondaryIdentifier", "symbol", "name", "sgdAlias", "phenotypeSummary",
    "chromosome.primaryIdentifier", "chromosomeLocation.start",
    "chromosomeLocation.end", "chromosomeLocation.strand"
    ]  
#"upstreamIntergenicRegion.chromosomeLocation.start","upstreamIntergenicRegion.chromosomeLocation.end","upstreamIntergenicRegion.length"
today = str(date.today())
filename = 'Dub_ORF_info_' + today +'.csv'
f = open(filename, 'w')

#,UpstreamIntergenicStart,UpstreamIntergenicEnd,UpstreamIntergenicLength
header_string = 'SystematicName,GeneName,FullName,SGD_Alias,Phenotype,ChrNum,Genome_Start,Genome_End,Strand\n'
# count = 1
# for view in views_list:    
#     if count == len(views_list):
#         header_string += view + '\n'
#     else:
#         header_string += view + ','
#     count += 1
#     print(count)
# print(header_string + 'new line')
          
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

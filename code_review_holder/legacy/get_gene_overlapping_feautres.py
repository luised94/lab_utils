# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:50:36 2020

@author: liusm
"""


from intermine.webservice import Service

yeast_url = "https://yeastmine.yeastgenome.org/yeastmine/service"
service = Service(yeast_url)

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

for row in query.rows():
    print(row["secondaryIdentifier"], row["symbol"], row["proteins.sequence.residues"], \
        row["proteins.sequence.length"], row["proteins.symbol"], row["locations.start"], \
        row["locations.end"], row["chromosome.primaryIdentifier"])





query2 = service.new_query("ARS")
query2.views
query2.add_view(
    "primaryIdentifier", "secondaryIdentifier", "symbol", "name", "locations.start", "locations.end", "chromosome.primaryIdentifier"
)
query2.add_constraint("ARS", "IN", "ARSs", code = "A")

for row in query2.rows():
    print(row["primaryIdentifier"], row["secondaryIdentifier"], row["symbol"], row["name"], \
          row["locations.start"], row["locations.end"])
        


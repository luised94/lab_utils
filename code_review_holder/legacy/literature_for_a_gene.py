#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "publicationAnnotations.publication.pubMedId",
    "publicationAnnotations.publication.citation",
    "publicationAnnotations.publication.year",
    "publicationAnnotations.literatureTopics.name", "primaryIdentifier",
    "symbol", "secondaryIdentifier", "sgdAlias", "name"
)

# This query's custom sort order is specified below:
query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "B")
query.add_constraint("Gene", "LOOKUP", "HEM2", code = "A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B")

for row in query.rows():
    print(row["publicationAnnotations.publication.pubMedId"], \
        row["publicationAnnotations.publication.citation"], \
        row["publicationAnnotations.publication.year"], \
        row["publicationAnnotations.literatureTopics.name"], row["primaryIdentifier"], \
        row["symbol"], row["secondaryIdentifier"], row["sgdAlias"], row["name"])


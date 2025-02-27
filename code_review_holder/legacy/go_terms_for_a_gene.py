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
    "primaryIdentifier", "secondaryIdentifier", "symbol", "featureType",
    "qualifier", "goAnnotation.ontologyTerm.identifier",
    "goAnnotation.ontologyTerm.name", "goAnnotation.ontologyTerm.namespace",
    "goAnnotation.evidence.code.code", "goAnnotation.qualifier",
    "goAnnotation.evidence.code.withText", "goAnnotation.annotationExtension",
    "goAnnotation.evidence.code.annotType",
    "goAnnotation.evidence.publications.pubMedId",
    "goAnnotation.evidence.publications.citation"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "F")
query.add_constraint("status", "IS NULL", code = "C")
query.add_constraint("status", "=", "Active", code = "B")
query.add_constraint("Gene", "LOOKUP", "YAL018C", code = "A")

# Your custom constraint logic is specified with the code below:
query.set_logic("(B or C) and F and A")

for row in query.rows():
    print(row["primaryIdentifier"], row["secondaryIdentifier"], row["symbol"], row["featureType"], \
        row["qualifier"], row["goAnnotation.ontologyTerm.identifier"], \
        row["goAnnotation.ontologyTerm.name"], row["goAnnotation.ontologyTerm.namespace"], \
        row["goAnnotation.evidence.code.code"], row["goAnnotation.qualifier"], \
        row["goAnnotation.evidence.code.withText"], row["goAnnotation.annotationExtension"], \
        row["goAnnotation.evidence.code.annotType"], \
        row["goAnnotation.evidence.publications.pubMedId"], \
        row["goAnnotation.evidence.publications.citation"])


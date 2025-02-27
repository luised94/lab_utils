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

# Type constraints should come early - before all mentions of the paths they constrain
query.add_constraint("interactions.participant2", "Gene")

# The view specifies the output columns
query.add_view(
    "primaryIdentifier", "symbol", "secondaryIdentifier", "sgdAlias", "name",
    "organism.shortName", "interactions.details.annotationType",
    "interactions.details.phenotype", "interactions.details.role1",
    "interactions.participant2.symbol",
    "interactions.participant2.secondaryIdentifier",
    "interactions.details.experiment.interactionDetectionMethods.identifier",
    "interactions.details.experiment.name",
    "interactions.details.relationshipType", "interactions.details.note"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "B")
query.add_constraint("Gene", "LOOKUP", "act1", code = "A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B")

for row in query.rows():
    print(row["primaryIdentifier"], row["symbol"], row["secondaryIdentifier"], row["sgdAlias"], \
        row["name"], row["organism.shortName"], row["interactions.details.annotationType"], \
        row["interactions.details.phenotype"], row["interactions.details.role1"], \
        row["interactions.participant2.symbol"], \
        row["interactions.participant2.secondaryIdentifier"], \
        row["interactions.details.experiment.interactionDetectionMethods.identifier"], \
        row["interactions.details.experiment.name"], row["interactions.details.relationshipType"], \
        row["interactions.details.note"])


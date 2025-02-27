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
query = service.new_query("Publication")

# The view specifies the output columns
query.add_view("pubMedId", "citation", "title", "journal", "literatureTopics.name")

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Publication.pubMedId", "ASC")

# You can edit the constraint values below
query.add_constraint("literatureTopics.name", "=", "Primary Literature", code = "A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A")

for row in query.rows():
    print(row["pubMedId"], row["citation"], row["title"], row["journal"], \
        row["literatureTopics.name"])


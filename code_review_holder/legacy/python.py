#!/usr/bin/env python2

import urllib2
import json

# Define request
acceptHeader = 'application/json' # text/csv and text/plain supported
request = urllib2.Request("https://mobidb.org/ws/P04050/uniprot", headers={"Accept" : acceptHeader})

# Send request
response = urllib2.urlopen(request)

# Parse JSON response di Python dict
data = json.load(response)

# handle data
print data


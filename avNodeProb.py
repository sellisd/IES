#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from pyies.userOptions import basePath
from collections import Counter
import os.path

# read spEvents into dictionary node->sp event
spEventsF = os.path.join(basePath, 'analysis', 'tables', 'spEvents1.dat')
avNodeProbF = os.path.join(basePath, 'analysis', 'tables', 'avNodeProb1.dat')
nodeDictF = os.path.join(basePath, 'analysis', 'tables', 'nodeDictionary1.dat')
homIESLinkF = os.path.join(basePath, 'analysis', 'tables', 'homIES1.columns.link')
# load node dictionary
rb2phyldog = {}
with open(nodeDictF, 'r') as f:
    f.readline()
    for line in f:
        line = line.rstrip()
        (cluster, r, phyldog, rb) = line.split("\t")
        rb2phyldog[(cluster, rb)] = phyldog

# load node id speciation events (ND:S)
node2Event = {}
with open(spEventsF, 'r') as f:
    f.readline()
    for line in f:
        line = line.rstrip()
        (cluster, nodeP, spEvent) = line.split("\t")
        node2Event[nodeP] = spEvent

# load column Ids homIES id
column2homIESid = {}
with open(homIESLinkF, 'r') as f:
    f.readline()
    for line in f:
        line = line.rstrip()
        (geneFamily, homIES, column) = line.split("\t")
        column2homIESid[(geneFamily, column)] = homIES

# read avnodProb use dictionary to translate from node to sp event
# sum for each sp event
avProb = Counter()
with open(avNodeProbF, 'r') as f:
    f.readline()
    for line in f:
        line = line.rstrip()
        (cluster, rb, iesColumn, presence) = line.split("\t")
        #homIESid = column2homIESid[(cluster, iesColumn)]
        nodeP = rb2phyldog[(cluster, rb)]
        print("\t".join([cluster, nodeP]))
        spEvent = node2Event[nodeP]
        avProb[spEvent] += float(presence)

print(avProb)

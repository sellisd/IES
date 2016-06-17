#!/usr/bin/python
from __future__ import print_function
from pyies.functions import *
""" In each homologous IES column of all gene families find the speciation nodes in which an IES was present with probability larger than cutoff"""

# find correspondance of phyldog nodedes to speciation events
speF = open('/home/dsellis/data/IES/analysis/tables/spEvents.dat', 'r')
v = Vividict() # node-sp.event correspondence
speF.readline()
for line in speF:
    (geneFamily, nodeP, spEvent) = line.split()
    v[geneFamily][nodeP] = spEvent

# read node dictionary and each phyldog tree and create list of nodeRb and speciation events
ndF = open('/home/dsellis/data/IES/analysis/tables/nodeDictionary.dat', 'r')
d = Vividict(); # all node dictionary
ndF.readline() # strip header
for line in ndF:
    line = line.rstrip()
    (geneFamily, r, phyldog, rb) = line.split()
    if v[geneFamily][phyldog]: # if a speciation node
        d[geneFamily][rb] = v[geneFamily][phyldog]

# read average node probabilities file
cutoff = 0.99
asrF = open('/home/dsellis/data/IES/analysis/tables/avNodeProb.dat', 'r')
header = asrF.readline()

# find all speciation nodes in which an IES was present with prob > cutoff
homies = {}
for line in asrF:
    line = line.rstrip()
    (geneFamily, nodeRb, iesColumn, presence) = line.split() # nodes are in rb format
    # find what type of speciation event nodeRb corresponds to
    if float(presence) > cutoff:
        spe = d[geneFamily][nodeRb]
        if spe:
            print(geneFamily, iesColumn, spe)


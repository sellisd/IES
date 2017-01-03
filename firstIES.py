#!/usr/bin/python
from __future__ import print_function
from pyies.functions import *
from ete3 import Tree
import re, sys
import os.path
from userOptions import basePath
# Find when an IES was inserted.

# In each homologous IES column of all gene families find the speciation nodes
# in which an IES was present with probability larger than cutoff.

asrRun = sys.argv[1]
# find correspondance of phyldog nodes to speciation events
speF = open(os.path.join(basePath, 'analysis', 'tables', 'spEvents'+ asrRun +'.dat'), 'r')
v = Vividict() # node-sp.event correspondence
speF.readline()
for line in speF:
    (geneFamily, nodeP, spEvent) = line.split()
    v[geneFamily][nodeP] = spEvent

# read node dictionary and each phyldog tree and create list of nodeRb and speciation events
ndF = open(os.path.join(basePath, 'analysis', 'tables', 'nodeDictionary' + asrRun + '.dat'), 'r')
d = Vividict(); # all node dictionary
ndF.readline() # strip header
for line in ndF:
    line = line.rstrip()
    (geneFamily, r, phyldog, rb) = line.split()
    if v[geneFamily][phyldog]: # if a speciation node
        d[geneFamily][rb] = v[geneFamily][phyldog]

# read average node probabilities file
cutoff = 0.99
asrF = open(os.path.join(basePath, 'analysis/', 'tables', 'avNodeProb' + asrRun + '.dat'), 'r')
header = asrF.readline()

# find all speciation nodes in which an IES was present with prob > cutoff
print('geneFamily' + "\t" + 'iesColumn' + "\t" + 'spEvent')
homies = {}
for line in asrF:
    line = line.rstrip()
    (geneFamily, nodeRb, iesColumn, presence) = line.split() # nodes are in rb format
    # find what type of speciation event nodeRb corresponds to
    if float(presence) > cutoff:
        spe = d[geneFamily][nodeRb]
        if spe:
            homies.setdefault((geneFamily,iesColumn), []).append(spe)

# load species tree
t = Tree(os.path.join(basePath, 'analysis/', 'phyldogT' + asrRun, 'results', 'OutputSpeciesTree_ConsensusNumbered.tree'), format = 1)

# replace species names with speciation events and add 0 to root node label
for l in t.traverse():
    if(l.is_leaf()):
        l.name = (re.sub(r'.+_.+_(\d+)', r'\1', l.name))
    elif(l.is_root()):
        if(l.name):
            quit(l.name) # is root named?
        else:
            l.name = '0'


for (geneFamily, iesColumn) in homies:
    L = homies[(geneFamily, iesColumn)]
    if(len(L) == 1):
        # if only one node
        print('\t'.join([geneFamily, iesColumn, L[0]]))
    else:
        ancestor = t.get_common_ancestor(L)
        print('\t'.join([geneFamily, iesColumn, ancestor.name]))

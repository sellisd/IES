#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree
import re
"""
Calculate all speciation node pairs on a PHYLDOG tree
"""

# Functions
def path2anc(nodeO):
    """
    Find path from root to node
    PARAMETERS: reference to a node object
    RETURN    : a list of node ids on the path from the root to the input node
    """
    path = []
    while nodeO.up:
        path.append(nodeO)
        nodeO = nodeO.up
    path.append(nodeO) # add root
    path.reverse()
    pathNames = []
    for i in path:
        pathNames.append(int(i.name))
    return pathNames

def nodePairs(l):
    """
    From a path of nodes generate all pairs
    """
    pairs = {}
    for i, itemI in enumerate(l):
        for itemJ in l[(i+1):]:
            pairs[(itemI,itemJ)] = 1
    return pairs

# load species tree
t = Tree('/home/dsellis/data/IES/analysis/phyldog/results/OutputSpeciesTree_ConsensusNumbered.tree', format = 1)


# replace species names with speciation events and add 0 to root node label
for l in t.traverse():
    if(l.is_leaf()):
        l.name = (re.sub(r'.+_.+_(\d+)', r'\1', l.name))
    elif(l.is_root()):
        if(l.name):
            quit(l.name) # is root named?
        else:
            l.name = '0'

allPairs = {}
#print(t.get_ascii(show_internal=True))
# get all paths to root
for l in t.get_leaves():
    p = path2anc(l)
    pairsInPath = nodePairs(p)
    allPairs.update(pairsInPath)

# sort
lp = sorted(allPairs.keys(), key=lambda tup: (tup[0], tup[1]))

# and print
for i in lp:
    print(i[0], i[1])



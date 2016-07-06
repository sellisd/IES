#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree
import re
""" for each IES find most recent common ancestor in species tree """

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


# read line by line firstIES.dat
f = open('/home/dsellis/data/IES/analysis/tables/firstIES.dat', 'r')
header = f.readline()
print('\t'.join(['cluster', 'iesColumn', 'spEvent']))
for line in f:
    line = line.rstrip()
    L = line.split()
    cluster = L.pop(0)
    iesColumn = L.pop(0)
    if(len(L) == 1):
        # if only one node
        print('\t'.join([cluster, iesColumn, L[0]]))
    else:
        ancestor = t.get_common_ancestor(L)
        print('\t'.join([cluster, iesColumn, ancestor.name]))

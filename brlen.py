#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from ete3 import Tree
import re

phyldogTreeFile = "/home/dsellis/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree"
brlenTreeFile = "/home/dsellis/data/IES/analysis/sgf/topoConstrSimple.treefile"
b = Tree(brlenTreeFile)
b.set_outgroup(b&"Tetrahymena_thermophila")
brlenD = {}
for node in b.traverse():
    leaveNames = [x.name for x in node.get_leaves()]
    leaveNames.sort()
    brlenD[tuple(leaveNames)] = x.dist

t = Tree(phyldogTreeFile)

def numbered2name(string):
    """Extract species name from PHYLDOG numbered leaf."""
    return(re.sub(r'(.+_.+)_\d+', r'\1', string))

nodeKey = {}
for node in t.traverse():
    PHYLDOGid = ''
    if node.is_leaf():
        PHYLDOGid = (re.sub(r'.+_.+_(\d+)', r'\1', node.name))
    elif node.is_root():
        PHYLDOGid = '0'
    else:
        PHYLDOGid = str(int(node.support))
    node.add_feature("PHYLDOGid", PHYLDOGid)
    leaveNames = [numbered2name(x.name) for x in node.get_leaves()]
    leaveNames.sort()
#    print(tuple(leaveNames))
    print(PHYLDOGid + ' ' + str(brlenD[tuple(leaveNames)]))
    node.dist = brlenD[tuple(leaveNames)]

t.show()
#Read a PHYLDOG output tree and a newick tree with branch lengths and merge information, print branches with branch lengths
#read phyldog tree, loop through nodes, for each node find all descendant leafs
##read branch length tree. loop through nodes, for each node find all descendant leafs
#find node name with same leaf ends and its branch length

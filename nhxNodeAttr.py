#!/usr/bin/python
from __future__ import print_function
from os import listdir
import os.path
from sys import argv
from ete3 import Tree

    """Parse node attributes from PHYLDOG reconciled tree files."""

baseP = argv[1]
treeF = [i for i in listdir(baseP) if i[-15:] == '.ReconciledTree']
for i in treeF:
    geneFamily = i[0:-15]
    t = Tree(os.path.join(baseP, i))
    for n in t.traverse():
        print(geneFamily, n.Ev, n.S, n.ND)

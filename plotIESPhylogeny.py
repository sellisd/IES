#!/usr/bin/python
from __future__ import print_function
import re
from ete3 import Tree, NodeStyle, TreeStyle
from pyies.functions import readPalette

t = Tree('/home/dsellis/data/IES/analysis/mies/tempfiles/cl35886.aln.contree')
cp = readPalette()

#cp  = ('1','2','3','4','5','6','7','8','9','10')
for node in t:
    m = re.search('(ppr|pbi|pte|ppe|pse|poc|ptr|pso|pca|tth)', node.name)
    if(m):
        nstyle = NodeStyle()
        nstyle["size"] = 10
        abr = m.group(0)
        if abr == 'ppr':
            nstyle['bgcolor'] = cp[0]
        elif abr == 'pbi':
            nstyle['bgcolor'] = cp[1]
        elif abr == 'pte':
            nstyle['bgcolor'] = cp[2]
        elif abr == 'ppe':
            nstyle['bgcolor'] = cp[3]
        elif abr == 'pse':
            nstyle['bgcolor'] = cp[4]
        elif abr == 'poc':
            nstyle['bgcolor'] = cp[5]
        elif abr == 'ptr':
            nstyle['bgcolor'] = cp[6]
        elif abr == 'pso':
            nstyle['bgcolor'] = cp[7]
        elif abr == 'pca':
            nstyle['bgcolor'] = cp[8]
        elif abr == 'tth':
            nstyle['bgcolor'] = cp[9]
        else:
            nstyle['bgcolor'] = "#000000"
        node.set_style(nstyle)


ts = TreeStyle()
#ts.show_leaf_name = False
ts.mode = 'c'
t.show(tree_style = ts)

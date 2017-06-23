#!/usr/bin/python
#figure IES loss
from __future__ import print_function
from ete3 import Tree, NodeStyle, TreeStyle
from pyies.userOptions import basePath
import os.path
import string

geneFamily = '63'
analysis = '2'
phyldogResultsPath = "analysis/phyldogT" + analysis + "/results"
treeF = os.path.join(basePath, 'analysis', 'phyldogT' + analysis, 'results', geneFamily + ".ReconciledTree"
)
t = Tree(treeF)
ts = TreeStyle()
ts.complete_branch_lines_when_necessary = False

nstyle = NodeStyle()
nstyle["size"] = 0
for n in t.traverse():
   n.set_style(nstyle)

t.show(tree_style = ts)
t.render("/Users/dsellis/projects/IES/src/figureData/tree." + geneFamily + ".svg", tree_style = ts)

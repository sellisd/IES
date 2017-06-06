#!/usr/bin/python
from ete3 import Tree, NodeStyle
from pyies.userOptions import basePath
import os.path

spTree1F = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree1.nhx')
spTree2F = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree2.nhx')
spTree3F = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree3b.nhx')
outT1 = os.path.join(basePath, 'analysis', 'figures', 'spTree1.png')
outT2 = os.path.join(basePath, 'analysis', 'figures', 'spTree2.png')
outT3 = os.path.join(basePath, 'analysis', 'figures', 'spTree3.png')

for tF, oF in zip([spTree1F, spTree2F, spTree3F], [outT1, outT2, outT3]):
    t = Tree(tF)
    t.sort_descendants(attr='name')
    for node in t.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        node.set_style(nstyle)
    #t.show()
    t.render(oF)

#!/usr/bin/evn python
from __future__ import print_function
from __future__ import division
from pyies.userOptions import basePath
from ete3 import TextFace, TreeStyle, NodeStyle
import pandas as pd
# only decorate a tree that has already PHYDLOG ids
outputF = ""
inputGLF = os.path.join(basePath, 'analysis', 'tables', 'gainLossSum1.dat')
inputTF = "/Volumes/WDC/data/IES/analysis/sgf/topoConstrSimple.treefile"
gl = pd.read_csv(inputGLF)
t = Tree(inputTF)

ts = TreeStyle()

# calculate branch colors
gainL = [] # list with all rates of gain
lossL = [] # list with all rates of loss

for node in t.iter_descendants():
    gainL.append(node.gain)
    lossL.append(node.loss)

bcrg = scaleCol(gl.pgain.tolist())  # Branch Colors for Rates of Gain
bcrl = scaleCol(gl.ploss.tolist())  # Branch Colors for Rates of Loss

# make a "gain" and a "loss" copy of the tree
tg = t.copy()
tl = t.copy()

for node in tg.iter_descendants(): # do not include root
    style = NodeStyle()
    gainString = "+%.2f" % (0.001*node.gain)
    style["vt_line_color"] = bcrg[node.gain]
    style["hz_line_color"] = bcrg[node.gain]
    style["size"] = 0
    node.add_face(TextFace(gainString), column = 0, position = "branch-top")
    node.set_style(style)

for node in tl.iter_descendants():
    style = NodeStyle()
    lossString = "-%.2f" % (10000*node.loss)
    style["vt_line_color"] = bcrl[node.loss]
    style["hz_line_color"] = bcrl[node.loss]
    style["size"] = 0
#    node.add_face(TextFace(node.PHYLDOGid), column = 0, position = "float")
    node.add_face(TextFace(lossString), column = 0, position = "branch-bottom")
    node.set_style(style)

if outputFileBaseName:
    tg.write(features = ["PHYLDOGid", "name", "gain"], outfile = outputFileBaseName + ".gain.tre")
    tl.write(features = ["PHYLDOGid", "name", "loss"], outfile = outputFileBaseName + ".loss.tre")
    if doNotDraw:
        print(tg.get_ascii(show_internal = True, attributes = ["PHYLDOGid", "name", "loss"]))
        print(tl.get_ascii(show_internal = True, attributes = ["PHYLDOGid", "name", "gain"]))
    else:
        tg.render(outputFileBaseName + ".gain.png", tree_style = ts)
        tl.render(outputFileBaseName + ".loss.png", tree_style = ts)
else:
    t.show()

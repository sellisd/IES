#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from pyies.userOptions import basePath
from pyies.functions import scaleCol
from ete3 import Tree, TextFace, TreeStyle, NodeStyle
import pandas as pd
import os.path

# only decorate a tree that has already PHYDLOG ids
for asrRun in (1,2,3):
    inputGLF = os.path.join(basePath, 'analysis', 'tables', 'gainLossSum' + str(asrRun) + '.dat')
    spTreeF  = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree' + str(asrRun) + '.nhx')
    outP     = os.path.join(basePath, 'analysis', 'figures', 'spTree' + str(asrRun))

    gl = pd.read_csv(inputGLF, sep = "\t")

    t = Tree(spTreeF)

    ts = TreeStyle()

    # calculate branch colors
    gainL = [] # list with all rates of gain
    lossL = [] # list with all rates of loss

    bcrg = scaleCol(gl.pgain.tolist())  # Branch Colors for Rates of Gain
    bcrl = scaleCol(gl.ploss.tolist())  # Branch Colors for Rates of Loss
    print(bcrg)
    quit()
    # make a "gain" and a "loss" copy of the tree
    tg = t.copy()
    tl = t.copy()
    for node in tg.iter_descendants(): # do not include root
        if node.up.is_root():
            pgain = gl.pgain[(gl.fromNode == 0) & (gl.toNode == int(node.ND))].tolist()
        else:
            pgain = gl.pgain[(gl.fromNode == int(node.up.ND)) & (gl.toNode == int(node.ND))].tolist()
        style = NodeStyle()
        pgain = pgain[0]
        gainString = "+%.2f" % (pgain)
        style["vt_line_color"] = bcrg[pgain]
        style["hz_line_color"] = bcrg[pgain]
        style["hz_line_width"] = 3
        style["vt_line_width"] = 3
        style["size"] = 0
        node.add_face(TextFace(gainString), column = 0, position = "branch-top")
        node.set_style(style)

    for node in tl.iter_descendants():
        if node.up.is_root():
            ploss = gl.ploss[(gl.fromNode == 0) & (gl.toNode == int(node.ND))].tolist()
        else:
            ploss = gl.ploss[(gl.fromNode == int(node.up.ND)) & (gl.toNode == int(node.ND))].tolist()
        ploss = ploss[0]
        style = NodeStyle()
        lossString = "-%.2f" % (ploss*100)
        style["vt_line_color"] = bcrl[ploss]
        style["hz_line_color"] = bcrl[ploss]
        style["hz_line_width"] = 3
        style["vt_line_width"] = 3
        style["size"] = 0
        node.add_face(TextFace(lossString), column = 0, position = "branch-bottom")
        node.set_style(style)
    tg.show()
    tg.render(outP + ".gain.png", tree_style = ts)
    tl.render(outP + ".loss.png", tree_style = ts)

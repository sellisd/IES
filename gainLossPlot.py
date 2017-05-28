#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from pyies.userOptions import basePath
from pyies.functions import scaleCol
from ete3 import Tree, TextFace, TreeStyle, NodeStyle
import pandas as pd
import os.path
import matplotlib.cm as cm
import matplotlib.colors as colors

# only decorate a tree that has already PHYDLOG ids
for asrRun in ('1','2','3'):
    inputGLF = os.path.join(basePath, 'analysis', 'tables', 'gainLossNormBrLen' + str(asrRun) + '.dat')
    outP     = os.path.join(basePath, 'analysis', 'figures', 'spTree' + str(asrRun))
    if asrRun == '3':
        analysisNumber = '3b'
    else:
        analysisNumber = asrRun
    spTreeF  = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree' + analysisNumber + '.nhx')
    gl = pd.read_csv(inputGLF, sep = "\t")
    t = Tree(spTreeF)
    ts = TreeStyle()
    # calculate branch colors
    gainL = [] # list with all rates of gain
    lossL = [] # list with all rates of loss
    gm = gl.rgain.min()
    gM = gl.rgain.max()
    lm = gl.rloss.min()
    lM = gl.rloss.max()
    #bcrg = scaleCol(gl.pgain.tolist())  # Branch Colors for Rates of Gain
    #bcrl = scaleCol(gl.ploss.tolist())  # Branch Colors for Rates of Loss
    # make a "gain" and a "loss" copy of the tree
    tg = t.copy()
    tl = t.copy()
    gcm = cm.ScalarMappable(norm = colors.Normalize(vmin = gm, vmax = gM), cmap = "coolwarm")
    lcm = cm.ScalarMappable(norm = colors.Normalize(vmin = lm, vmax = lM), cmap = "coolwarm")
    for node in tg.iter_descendants(): # do not include root
        if node.up.is_root():
            rgain = gl.loc[(gl.fromNode == 0) & (gl.toNode == int(node.ND)), 'rgain']
        else:
            rgain = gl.loc[(gl.fromNode == int(node.up.ND)) & (gl.toNode == int(node.ND)), 'rgain']
        if rgain.empty:
            continue
        rgain = rgain.item()
        style = NodeStyle()
        gainString = "+%.2f" % (rgain)
        #pick colors
        ci = colors.rgb2hex(gcm.to_rgba(rgain)[:3])
        style["vt_line_color"] = ci
        style["hz_line_color"] = ci
        style["hz_line_width"] = 3
        style["vt_line_width"] = 3
        style["size"] = 0
        node.add_face(TextFace(gainString), column = 0, position = "branch-top")
        node.set_style(style)

    for node in tl.iter_descendants():
        if node.up.is_root():
            rloss = gl.rloss[(gl.fromNode == 0) & (gl.toNode == int(node.ND))]
        else:
            rloss = gl.rloss[(gl.fromNode == int(node.up.ND)) & (gl.toNode == int(node.ND))]
        if rloss.empty:
            continue
        rloss = rloss.item()
        style = NodeStyle()
        lossString = "-%.2f" % (rloss)
        ci = colors.rgb2hex(lcm.to_rgba(rloss)[:3])
        style["vt_line_color"] = ci
        style["hz_line_color"] = ci
        style["hz_line_width"] = 3
        style["vt_line_width"] = 3
        style["size"] = 0
        node.add_face(TextFace(lossString), column = 0, position = "branch-bottom")
        node.set_style(style)
    tg.show()
    tl.show()
    tg.render(outP + ".gain.png", tree_style = ts)
    tl.render(outP + ".loss.png", tree_style = ts)

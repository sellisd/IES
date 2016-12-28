#!/usr/bin/python
from __future__ import division
from __future__ import print_function
from ete3 import TextFace, NodeStyle
from pyies.functions import phyldogSpeciesTree
import sys, getopt
import os.path
import matplotlib
import matplotlib.cm
#/phyldogEventSum.py pathToPhyldog/results/*_Events.txt > ~/data/IES/analysis/phyldogT1/geneDuplLossT2.dat

def scaleCol(L):
    """ Create color scale from list."""
    m = matplotlib.cm.ScalarMappable(cmap = "coolwarm")
    R = {}
    for (i, rgba) in enumerate(m.to_rgba(L)):
        rgb = rgba[:3]
        R[L[i]] = matplotlib.colors.rgb2hex(rgb)
    return(R)

analysis = '2'
basePath = "/Users/dsellis/data/IES/"
outputFile = ""
normalize = 0
colorBy = "d"

usage = """
usage:

./phyldogEventPlot.py [OPTIONS]

    where OPTIONS can be any of the following:
    -b <basePath> Default /Users/dsellis/data/IES
    -a [1|2|3] Choice of species tree analysis to use (Default: 2)
    -o outputFile to render image, if omitted visualize in GUI
    -c [d|l] color branches by duplications or losses (number of events or rate is determined by  -n)
    -n if provided normalize by branch length, else print numbers
    -h this help screen

The coding for species tree analysis is:
1. PHYLDOG species tree
2. concatenate with simple evolution model
3. concatenate with best codon models

Example:
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -c d -a 1 -o ~/data/IES/analysis/figures/geneConvRates1.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -c l -a 1 -o ~/data/IES/analysis/figures/geneLossRates1.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -c d -a 2 -o ~/data/IES/analysis/figures/geneConvRates2.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -c l -a 2 -o ~/data/IES/analysis/figures/geneLossRates2.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -c d -a 3 -o ~/data/IES/analysis/figures/geneConvRates3.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -c l -a 3 -o ~/data/IES/analysis/figures/geneLossRates3.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py    -c d -a 1 -o ~/data/IES/analysis/figures/geneConvEvents1.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py    -c l -a 1 -o ~/data/IES/analysis/figures/geneLossEvents1.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py    -c d -a 2 -o ~/data/IES/analysis/figures/geneConvEvents2.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py    -c l -a 2 -o ~/data/IES/analysis/figures/geneLossEvents2.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py    -c d -a 3 -o ~/data/IES/analysis/figures/geneConvEvents3.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py    -c l -a 3 -o ~/data/IES/analysis/figures/geneLossEvents3.png &
"""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hb:a:o:c:n")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == "-b":
        basePath = arg
    elif opt == "-a":
        if arg in ['1','2','3']:
            analysis = arg
        else:
            print("Unknown value for option analysis (-a)")
            print(usage)
            sys.exit(1)
    elif opt == "-o":
        outputFile = arg
    elif opt == "-n":
        normalize = 1
    elif opt == "-c":
        if arg in ['d', 'l']:
            colorBy = arg
        else:
            print("Unknown value for option color by (-c)")
            print(usage)
            sys.exit(1)

brlenTreeFile = ""
if analysis == '1':
    brlenTreeFile  = os.path.join(basePath, "analysis/sgf/topoConstrSimple.treefile")
elif analysis == '2':
    brlenTreeFile = os.path.join(basePath, "analysis/phyldogT2/concatSimple.r.tre")
elif analysis == '3':
    brlenTreeFile = os.path.join(basePath, "analysis/phyldogT3/concat.r.tre")
else:
    sys.exit(1)

eventsfile = os.path.join(basePath, 'analysis/phyldogT1/geneDuplLossT'+ analysis +'.dat') # all geneDuplLossT files are in phyldogT1
phyldogTreeFile = os.path.join(basePath, 'analysis/phyldogT' + analysis + '/results/OutputSpeciesTree_ConsensusNumbered.tree')
outgroupName = "Tetrahymena_thermophila"

t = phyldogSpeciesTree(phyldogTreeFile, brlenTreeFile, outgroupName)

brattr = {}
f = open(eventsfile, 'r')
f.readline() #get rid of header
for line in f:
    line = line.rstrip()
    (speciesNode, duplications, losses) = line.split("\t")
    brattr[speciesNode] = (duplications, losses)

# add rates on branch attribute table
for node in t.traverse():
    if node.PHYLDOGid in brattr:
        dEvents = brattr[node.PHYLDOGid][0] # duplication or gene conversion events
        lEvents = brattr[node.PHYLDOGid][1]
        if node.dist > 10**-10:
            dRate = int(dEvents) / node.dist
            lRate = int(lEvents) / node.dist
        else:
            dRate = 0 # NA
            lRate = 0
        brattr[node.PHYLDOGid] = (dEvents, lEvents, dRate, lRate)

# calculate branch colors
DL = [] # list with all duplications
LL = [] # list with all losses
DRL = []
LRL = []

for k in brattr:
    (duplications, losses, dRate, lRate) = brattr[k]
    DL.append(duplications)
    LL.append(losses)
    DRL.append(dRate)
    LRL.append(lRate)

bcde = scaleCol(DL)  # Branch Colors for Duplication Events
bcle = scaleCol(LL)  # Branch Colors for Loss Events
bcdr = scaleCol(DRL) # Branch Colors for Duplication Rates
bclr = scaleCol(LRL) # Branch Colors for Loss Rates

for node in t.traverse():
    if node.PHYLDOGid in brattr:
        if normalize:
            dRates = brattr[node.PHYLDOGid][2] # duplication or gene conversion events
            lRates = brattr[node.PHYLDOGid][3]
            style = NodeStyle()
            if colorBy == "d":
                style["vt_line_color"] = bcdr[dRates]
                style["hz_line_color"] = bcdr[dRates]
                node.add_face(TextFace(dRates), column = 0, position = "branch-top")
            elif colorBy == "l":
                style["vt_line_color"] = bclr[lRates]
                style["hz_line_color"] = bclr[lRates]
                node.add_face(TextFace(lRates), column = 0, position = "branch-bottom")
            node.set_style(style)
        else:
            dEvents = brattr[node.PHYLDOGid][0] # duplication or gene conversion events
            lEvents = brattr[node.PHYLDOGid][1]
            style = NodeStyle()
            if colorBy == "d":
                style["vt_line_color"] = bcde[dEvents]
                style["hz_line_color"] = bcde[dEvents]
                node.add_face(TextFace(dEvents), column = 0, position = "branch-top")
            elif colorBy == "l":
                style["vt_line_color"] = bcle[lEvents]
                style["hz_line_color"] = bcle[lEvents]
                node.add_face(TextFace(lEvents), column = 0, position = "branch-bottom")
            node.set_style(style)

if outputFile:
    t.render(outputFile)
else:
    t.show()
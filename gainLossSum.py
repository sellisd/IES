#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from collections import Counter, defaultdict
from pyies.functions import phyldogSpeciesTree, scaleCol
import sys, getopt
from ete3 import TextFace, TreeStyle, NodeStyle
from decimal import *

# Normalize loss rate by total length of conserved blocks of alignments and both insertion and loss rate by branch lengths.

# program options
gainLossFile = "/home/dsellis/data/IES/analysis/tables/gainLoss1.dat"
gbFile = "/home/dsellis/data/IES/analysis/tables/gblocks.dat" # Gblocks file
brlenFile = "/home/dsellis/data/IES/analysis/sgf/topoConstrSimple.treefile"
phyldogTreeFile = "/home/dsellis/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree"
outgroupName = "Tetrahymena_thermophila"
outputFileBaseName = ""
doNotDraw = 0
normBrLen = 0

usage = """
usage:

./gainLossSum.py [OPTIONS]

where OPTIONS can be any of the following:

    -g: gainLossFile
    -b: gblocksFile
    -l: Newick File of tree with branch lengths
    -p: PHYLDOG outputFile
    -n: outgroup name (Default: Tetrahymena_thermophila)
    -o: output File Base Name (if not provided show in interactive tree viewer)
    -d: do not draw tree, output text file and print ASCII tree
    -r: normalize by branch lengths
    -h: this help screen
""";

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:b:l:p:n:o:d:r")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == '-g':
        gainLossFile = arg
    elif opt == '-b':
        gbFile = arg
    elif opt == '-l':
        brlenFile = arg
    elif opt == '-p':
        phyldogTreeFile = arg
    elif opt == '-n':
        outgroupName = arg
    elif opt == '-o':
        outputFileBaseName = arg
    elif opt == '-d':
        doNotDraw = 1
    elif opt == '-r':
        normBrLen = 1

# load gblock sizes and sum for each gene family
print("sum gblock sizes")
gf = open(gbFile, 'r')
gb = Counter()
for line in gf:
    line = line.rstrip()
    (geneFamily, begin, end) = line.split()
    gb[geneFamily] += int(end) - int(begin) + 1


# sum gain and loss probabilities along branches
print("sum gain and loss probability along paths")
gl = open(gainLossFile, 'r')
gl.readline() # header
sumgain = Counter()
nogain  = Counter()
sumloss = Counter()
noloss  = Counter()
Nij = defaultdict(set) # number of gene families with Si-Sj path

for line in gl:
    line = line.rstrip()
    (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) = line.split()
    sumloss[(fromNode, toNode)] += float(panc) * float(loss) # normalize rate of loss by the probability of being present
    noloss[(fromNode, toNode)] += 1
    sumgain[(geneFamily, fromNode, toNode)] += float(gain)
    nogain[(geneFamily, fromNode, toNode)] += 1
    Nij[(fromNode, toNode)].add(geneFamily)

# normalize rate of gain by gblocks length and number of paths
print("normalize")
pgain = Counter()
for k in sumgain:
    pgain[(k[1], k[2])] += sumgain[k] * gb[k[0]] / nogain[k]

# normalize rate of loss by number of paths
ploss = Counter()
for k in sumloss:
    ploss[k] = sumloss[k] / noloss[k]

# normalize rate of gain and loss by number of gene families with Si-Sj path
for k in pgain:
    pgain[k] /= len(Nij[k])
    ploss[k] /= len(Nij[k])


# normalize by branch length
t = phyldogSpeciesTree(phyldogTreeFile, brlenFile, outgroupName)
for k in pgain:
    node = t.search_nodes(PHYLDOGid=k[1])[0]
    if normBrLen == 1:
        pgain[k] /= float(node.dist)
    node.add_feature("gain", pgain[k])

for k in ploss:
    node = t.search_nodes(PHYLDOGid=k[1])[0]
    if normBrLen == 1:
        ploss[k] /= float(node.dist)
    node.add_feature("loss", ploss[k])

ts = TreeStyle()

# calculate branch colors
gainL = [] # list with all rates of gain
lossL = [] # list with all rates of loss

for node in t.iter_descendants():
    gainL.append(node.gain)
    lossL.append(node.loss)

bcrg = scaleCol(gainL)  # Branch Colors for Rates of Gain
bcrl = scaleCol(lossL)  # Branch Colors for Rates of Loss

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

#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from collections import Counter, defaultdict
from pyies.functions import phyldogSpeciesTree
import sys, getopt
from ete3 import TextFace, TreeStyle
from decimal import *

# Normalize loss rate by  total length of conserved blocks of alignments and both insertion and loss rate by branch lengths.

# program options
gainLossFile = "/home/dsellis/data/IES/analysis/tables/gainLoss1.dat"
gbFile = "/home/dsellis/data/IES/analysis/tables/gblocks.dat" # Gblocks file
brlenFile = "/home/dsellis/data/IES/analysis/sgf/topoConstrSimple.treefile"
phyldogTreeFile = "/home/dsellis/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree"
outgroupName = "Tetrahymena_thermophila"
outputFileFig = ""
outputFileNewick = ""

usage = "./gainLossSum.py -g <gainLossFile> -b <gblocksFile> -l <treeFile> -p <PHYLDOGoutputFile> -n <outgroupName> -o <outputFileFig> -w <outputFileNewick>"

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:b:l:p:n:o:w:")
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
        outputFileFig = arg
    elif opt == '-w':
        outputFileNewick = arg

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
    pgain[k] /= float(node.dist)
    node.add_feature("gain", pgain[k])

for k in ploss:
    node = t.search_nodes(PHYLDOGid=k[1])[0]
    ploss[k] /= float(node.dist)
    node.add_feature("loss", ploss[k])

ts = TreeStyle()

#print(t.get_ascii(show_internal = True, attributes = ["PHYLDOGid", "name", "loss"]))
#print(t.get_ascii(show_internal = True, attributes = ["PHYLDOGid", "name", "gain"]))

for node in t.iter_descendants(): # do not include root
    gainString = "+%d" % (node.gain)
    lossString = "-%.2f" % (10000*node.loss)
    node.add_face(TextFace(node.PHYLDOGid), column = 0, position = "float")
    node.add_face(TextFace(gainString), column = 0, position = "branch-top")
    node.add_face(TextFace(lossString), column = 0, position = "branch-bottom")

t.write(features = ["PHYLDOGid", "name", "loss"], outfile = outputFileNewick)
t.render(outputFileFig, tree_style = ts)
#t.show()

#!/usr/bin/python
from __future__ import print_function
from collections import Counter
from __future__ import division


# normalize loss rate by  total length of conserved blocks of alignments and both insertion and loss rate by branch lengths

gainLossFile = "/home/dsellis/data/IES/analysis/tables/gainLoss1.dat"
gbFile = "/home/dsellis/data/IES/analysis/tables/gblocks.dat" # Gblocks file
brlenFile = "/home/dsellis/data/IES/analysis/sgf/topoConstrSimple.treefile"
#/home/dsellis/data/IES/analysis/sgf/
# manually create table with PHYLDOG from to speciation tree nodes and corresponding branch length in each tree
# or make .key file for each OutputSpeciesTree_ConsensusNumbered
# and for brlenFiles
# and then merge to create table
# or read brlenFile
# and tree with species numbers
# for each node in species tree find corresponding in brlentree
# or seaview PHYLDOG with numbers
# and plot in R tree

# read brlenfile
# read PHYLDOG with numbers
# and make correspondences

#gainLossFile has in PHYLDOG node numbering
#brlenFile doesn't have any numbering
#TODO for each brlen branch find the corresponding PHYLDOG number
# write function that reads newick from PHYLDOG
# calculate total length of conserved blocks of alignments
gf = open(gbFile, 'r')
gb = Counter()
for line in gf:
    line = line.rstrip()
    (geneFamily, begin, end) = line.split()
    gb[geneFamily] += int(end) - int(begin) + 1


# load calculate branch lengths
bf = open(brlenFile, 'r')
bf.readline() # header
sumgain = Counter()
nogain  = Counter()
sumloss = Counter()
noloss  = Counter()
for line in bf:
    line = line.rstrip()
    (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) = line.split()
    sumloss[(fromNode, toNode)] += float(panc) * float(loss) # normalize rate of loss by the probability of being present
    noloss[(fromNode, toNode)] += 1
    # normalize rate of gain by total length of conserved blocks
    sumgain[(geneFamily, fromNode, toNode)] += float(gain)
    nogain[(geneFamily, fromNode, toNode)] += 1

pgain = Counter()
for k in sumgain:
    pgain[(k[1], k[2])] += sumgain[k] * nogain[k] /gb[k[0]]

ploss = Counter()
for k in sumloss:
    ploss[k] = sumloss[k] / noloss[k]

#normalize by branch length

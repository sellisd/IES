#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from collections import Counter
from pyies.functions import phyldogSpeciesTree

# normalize loss rate by  total length of conserved blocks of alignments and both insertion and loss rate by branch lengths

gainLossFile = "/home/dsellis/data/IES/analysis/tables/gainLoss1.dat"
gbFile = "/home/dsellis/data/IES/analysis/tables/gblocks.dat" # Gblocks file
brlenFile = "/home/dsellis/data/IES/analysis/sgf/topoConstrSimple.treefile"
phyldogTreeFile = "/home/dsellis/data/IES/analysis/phyldogT1/results/OutputSpeciesTree_ConsensusNumbered.tree"
outgroupName = "Tetrahymena_thermophila"

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
for line in gl:
    line = line.rstrip()
    (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) = line.split()
    sumloss[(fromNode, toNode)] += float(panc) * float(loss) # normalize rate of loss by the probability of being present
    noloss[(fromNode, toNode)] += 1
    sumgain[(geneFamily, fromNode, toNode)] += float(gain)
    nogain[(geneFamily, fromNode, toNode)] += 1

# normalize probability of gain by gblocks length
print("normalize")
pgain = Counter()
for k in sumgain:
    pgain[(k[1], k[2])] += sumgain[k] * nogain[k] / gb[k[0]]

ploss = Counter()
for k in sumloss:
    ploss[k] = sumloss[k] / noloss[k]

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

print(t.get_ascii(show_internal = True, attributes = ["PHYLDOGid", "name", "loss"]))
print(t.get_ascii(show_internal = True, attributes = ["PHYLDOGid", "name", "gain"]))

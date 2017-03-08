#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from collections import Counter
from pyies.userOptions import basePath
import os.path
import sys, getopt
from ete3 import Tree

# Normalize loss rate by total length of conserved blocks of alignments and both insertion and loss rate by branch lengths.

# program options
spNodePairsF       = ""
gainLossFile       = ""
gbFile             = ""
spTreeF            = ""
outgroupName       = ""
outputF            = ""
includeGF          = ""
nodePathsF         = ""
usage = """
usage:

./gainLossSum.py [OPTIONS]

where OPTIONS can be any of the following:

    -k: nodePaths file
    -g: gainLossFile
    -b: gblocksFile
    -t: species tree file with branch lengths
    -o: output file name
    -i: file with gene families to include in the analysis. Default all included
    -h: this help screen
""";

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:b:t:o:k:i")
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
    elif opt == '-t':
        spTreeF = arg
    elif opt == '-n':
        outgroupName = arg
    elif opt == '-o':
        outputF = arg
    elif opt == '-i':
        includeGF = arg
    elif opt == '-k':
        nodePathsF = arg

# if parameter defined load a list of gene families to include from the analysis
includedGeneFamilies = []
if includeGF:
    with open(includeGF, 'r') as f:
        for line in f:
            includedGeneFamilies = [line.rstrip() for line in f]

# load gblock sizes and sum for each gene family
print("sum gblock sizes")
gf = open(gbFile, 'r')
gb = Counter()
for line in gf:
    line = line.rstrip()
    (geneFamily, begin, end) = line.split()
    if (includeGF and (geneFamily in includedGeneFamilies)) or (not includeGF):
        gb[geneFamily] += int(end) - int(begin) + 1

# load kij
kij = Counter() # The number of paths connecting nodes i and j per gene family
with open(nodePathsF, 'r') as f:
    for line in f:
        line = line.rstrip()
        (cluster, fromNode, toNode, path) = line.split() # nodes in phyldog notation
        kij[(cluster, fromNode, toNode)] += 1

# sum gain and loss probabilities along branches
print("sum gain and loss probability along paths")
gl = open(gainLossFile, 'r')
gl.readline() # header
sumgain = Counter()
sumloss = Counter()
pgain = Counter()
ploss = Counter()
Ig = Counter()
iIg = Counter() # number of IES columns per gene family

for line in gl:
    line = line.rstrip()
    (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) = line.split()
    if (includeGF and (geneFamily in includedGeneFamilies)) or (not includeGF):
        sumgain[(geneFamily, fromNode, toNode)] += float(gain)
        sumloss[(fromNode, toNode)] += float(loss)
        iIg[(geneFamily, iesColumn)] = 1

for i in iIg:
    Ig[i[0]] +=1

# normalize by number of paths and for rate of gain also by gblocks length in nt
print("normalize")
for k in sumgain:
    pgain[(k[1], k[2])] += sumgain[k] / (gb[k[0]] * kij[k])
    ploss[(k[1], k[2])] += sumloss[k] / ( Ig[k[0]] * kij[k])

# normalize by branch length
# if a branch is not present in the species tree (e.g. skips a speciation events
# then in order to properly normalize we would have to add up the total branch
# as we don't need it for the figures we only keep the branch lengths observed
# in the species trees)
obspgain = {}
obsploss = {}
t = Tree(spTreeF)
for node in t.iter_descendants():
    fromNode = node.up.ND
    k = (fromNode, node.ND)
    obspgain[k]  = pgain[k]/float(node.dist)
    obsploss[k] = ploss[k]/float(node.dist)

with open(outputF, 'w') as fout:
    fout.write("\t".join(["fromNode", "toNode", "pgain", "ploss\n"]))

    for k in obspgain:
        fout.write("\t".join([str(k[0]), str(k[1]), str(obspgain[k]), str(obsploss[k]) + "\n"]))

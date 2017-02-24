#!/usr/bin/python
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
usage = """
usage:

./gainLossSum.py [OPTIONS]

where OPTIONS can be any of the following:

    -g: gainLossFile
    -b: gblocksFile
    -t: species tree file with branch lengths
    -o: output file name
    -i: file with gene families to include in the analysis. Default all included
    -h: this help screen
""";

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:b:t:o:i")
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

# sum gain and loss probabilities along branches
print("sum gain and loss probability along paths")
gl = open(gainLossFile, 'r')
gl.readline() # header
sumgain = Counter()
nogain  = Counter()
sumloss = Counter()
noloss  = Counter()
pgain = Counter()
ploss = Counter()

for line in gl:
    line = line.rstrip()
    (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) = line.split()
    if (includeGF and (geneFamily in includedGeneFamilies)) or (not includeGF):
        sumloss[(fromNode, toNode)] += float(panc) * float(loss) # normalize rate of loss by the probability of being present
        noloss[(fromNode, toNode)] += 1
        sumgain[(geneFamily, fromNode, toNode)] += float(gain)
        nogain[(fromNode, toNode)] += 1

# normalize rate of gain by gblocks length
print("normalize")
for k in sumgain:
    pgain[(k[1], k[2])] += sumgain[k] * gb[k[0]]

# normalize rate of gain and loss by total number of Si-Sj paths
for k in pgain:
    pgain[k] /= nogain[k]
    ploss[k] = sumloss[k] / noloss[k]

# normalize by branch length
t = Tree(spTreeF)
for k in pgain:
    node = t.search_nodes(ND=k[1])[0]
    pgain[k] /= float(node.dist)

for k in ploss:
    node = t.search_nodes(ND=k[1])[0]
    ploss[k] /= float(node.dist)

with open(outputF, 'w') as fout:
    fout.write("\t".join(["fromNode", "toNode", "pgain", "ploss\n"]))

    for k in pgain:
        fout.write("\t".join([k[0], k[1], str(pgain[k]), str(ploss[k]) + "\n"]))

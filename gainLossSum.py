#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import pandas as pd
from collections import Counter, defaultdict
from pyies.userOptions import basePath
import os.path
import sys, getopt
from ete3 import Tree

# Normalize loss rate by total length of conserved blocks of alignments and both insertion and loss rate by branch lengths.

# program options
spNodePairsF       = ""
gainLossF          = ""
gbFile             = ""
spTreeF            = ""
speciesNodePairsF  = ""
outgroupName       = ""
outputF            = ""
includeGF          = ""
nodePathsF         = ""
usage = """
usage:

./gainLossSum.py [OPTIONS]

where OPTIONS can be any of the following:

    -k: nodePaths file
    -g: gainLoss file
    -b: gblocksFile
    -t: species tree file with branch lengths
    -o: output file name
    -i: file with gene families to include in the analysis. Default all included
    -p: species node pairs file
    -h: this help screen
""";

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:b:t:o:k:p:i")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == '-g':
        gainLossF = arg
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
    elif opt == '-p':
        speciesNodePairsF = arg
    else:
        print("Unknown argument" + arg)
        print(usage)
        sys.exit(1)

npF = pd.read_table(nodePathsF, sep = "\t", index_col = False)
glF = pd.read_table(gainLossF, sep = "\t", index_col = False)
gbF = pd.read_table(gbFile, sep = "\t", index_col = False)
snpF = pd.read_table(speciesNodePairsF, sep = "\s+", index_col = False)
# if parameter defined load a list of gene families to include from the analysis
includedGeneFamilies = []
if includeGF:
    with open(includeGF, 'r') as f:
        for line in f:
            includedGeneFamilies = [line.rstrip() for line in f]

# load gblock sizes and sum for each gene family
print("sum gblock sizes")
gb = Counter() #ng
for (geneFamily, begin, end) in gbF.itertuples(index = False, name = None):
        gb[geneFamily] += int(end) - int(begin) + 1

kij = Counter() # The number of paths connecting nodes i and j per gene family
for (geneFamily, fromNode, toNode, path) in npF.itertuples(index = False, name = None):
    kij[(cluster, fromNode, toNode)] += 1 # nodes are in phyldog notation

print("sum gain and loss probability along paths")
sumgain = Counter() # pcij
Ig = defaultdict(set) # IES columns per gene family

df = pd.dataFrame(0, index =range(len(glF)), columns=["geneFamily", "pgij", "kij", "ng", "plij"])
#fill in data frame .
#  read by row gailloss.dat file and add to the correect row of df the rate of gain or loss
for (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) in glF.itertuples(index = False, name = None):
    df.loc[df["geneFamily"]==geneFamily,"pgij"] += pgain
    df.loc[df["geneFamily"]==geneFamily,"plij"] += ploss
    df.loc[df["geneFamily"]==geneFamily,"kij"]  = kij[(geneFamily, fromNode, toNode)]
    df.loc[df["geneFamily"]==geneFamily,"ng"]   = gb[geneFamily]
    Ig[geneFamily].add(iesColumn)

#print("\t".join(["geneFamily", "pcij", "kij", "ng"]))

#for k in sumgain:
#    print("\t".join([str(geneFamily), str(sumgain[k]), str(kij[k]), str(gb[k[0]])]))

#for each pair of speciation nodes perform a sum

quit()
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

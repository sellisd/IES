#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import pandas as pd
from collections import Counter, defaultdict
from pyies.userOptions import basePath
import os.path
import sys, getopt

# Calculate gain and loss rate.

# program options
spNodePairsF       = ""
gainLossF          = ""
gbFile             = ""
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
    else:
        print("Unknown argument" + arg)
        print(usage)
        sys.exit(1)

npF = pd.read_table(nodePathsF, sep = "\t", index_col = False)
glF = pd.read_table(gainLossF, sep = "\t", index_col = False)
gbF = pd.read_table(gbFile, sep = "\t", index_col = False)
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
    kij[(geneFamily, fromNode, toNode)] += 1 # nodes are in phyldog notation

print("sum gain and loss probability along paths")
sumgain = Counter() # pcij
sumloss = Counter()
Ig = defaultdict(set) # IES columns per gene family

#fill a data frame is too slow
#df = pd.DataFrame(0, index =range(len(glF)), columns=["geneFamily", "pgij", "kij", "ng", "plij"])
#  read by row gailloss.dat file and add to the correect row of df the rate of gain or loss
for (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) in glF.itertuples(index = False, name = None):
    sumgain[(geneFamily, fromNode, toNode)] += gain
    sumloss[(geneFamily, fromNode, toNode)] += loss
    #df.loc[df["geneFamily"]==geneFamily,"pgij"] += gain
    #df.loc[df["geneFamily"]==geneFamily,"plij"] += loss
    #df.loc[df["geneFamily"]==geneFamily,"kij"]  = kij[(geneFamily, fromNode, toNode)]
    #df.loc[df["geneFamily"]==geneFamily,"ng"]   = gb[geneFamily]
    Ig[geneFamily].add(iesColumn)

with open(outputF, 'w') as f:
    f.write("\t".join(["geneFamily", "fromNode", "toNode", "pcijGain", "pcijLoss", "kij", "ng\n"]))
    for k in sumgain:
        f.write("\t".join([str(geneFamily), str(k[1]), str(k[2]), str(sumgain[k]), str(sumloss[k]), str(kij[k]), str(gb[k[0]]) + "\n"]))

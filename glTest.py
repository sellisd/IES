#!/usr/bin/env python
# gl test
from __future__ import print_function
from __future__ import division
import pandas as pd
import os.path
from pyies.userOptions import basePath
from collections import Counter, defaultdict

# ./glTest.py > /
npF = pd.read_table(os.path.join(basePath, "analysis/tables/nodePaths2.dat"), sep = "\t", index_col = False)
glF = pd.read_table(os.path.join(basePath, "analysis/tables/gainLoss2.dat"), sep = "\t", index_col = False)
gbF = pd.read_table(os.path.join(basePath, "analysis/tables/gblocks.dat"), sep = "\t", index_col = False)
#  -t /Volumes/WDC/data/IES/analysis/iesdb/speciesTree2.nhx
gfF = pd.read_table(os.path.join(basePath, "analysis/iesdb/geneFamilydb.dat"), sep = "\t", index_col = False)
# Save as output a table with columns: genefamily, pgij kgij ng
#
# make a list of gene families in our analysis

gb = Counter() #ng
for (geneFamily, begin, end) in gbF.itertuples(index = False, name = None):
        gb[geneFamily] += int(end) - int(begin) + 1

sumgain = Counter() # pcij
Ig = defaultdict(set) # IES columns per gene family

for (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) in glF.itertuples(index = False, name = None):
    sumgain[(geneFamily, fromNode, toNode)] += gain
    Ig[geneFamily].add(iesColumn)

kij = Counter() # The number of paths connecting nodes i and j per gene family
for (geneFamily, fromNode, toNode, path) in npF.itertuples(index = False, name = None):
    kij[(geneFamily, fromNode, toNode)] += 1 # nodes are in phyldog notation

print("\t".join(["geneFamily", "fromNode", "toNode", "pgij", "kij", "ng"]))

for k in sumgain:
    print("\t".join([str(k[0]), str(k[1]), str(k[2]),str(sumgain[k]), str(kij[k]), str(gb[k[0]])]))

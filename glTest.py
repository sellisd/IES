#!/usr/bin/env python
# gl test
import pandas as pd

npF = pd.read_table("/Volumes/WDC/data/IES/analysis/tables/nodePaths2.dat", sep = "\t", index_col = False)
glF = pd.read_table("/Volumes/WDC/data/IES/analysis/tables/gainLoss2.dat", sep = "\t", index_col = False)
gbF = pd.read_table("/Volumes/WDC/data/IES/analysis/tables/gblocks.dat", sep = "\t", index_col = False)
#  -t /Volumes/WDC/data/IES/analysis/iesdb/speciesTree2.nhx
gfF = pd.read_table("/Volumes/WDC/data/IES/analysis/iesdb/geneFamilydb.dat", sep = "\t", index_col = False)
# make a list of gene families in our analysis
gb = Counter()
for (geneFamily, begin, end) in gbF.itertuples(index = False, name = None):
        gb[geneFamily] += int(end) - int(begin) + 1

sumgain = Counter()
Ig = defaultdict(set) # IES columns per gene family

for (geneFamily, iesColumn, fromNode, toNode, panc, gain, loss) in gl.itertuples(index = False, name = None):
    sumgain[(geneFamily, fromNode, toNode)] += gain
    Ig[geneFamily].add(iesColumn)

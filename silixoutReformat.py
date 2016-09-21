#!/usr/bin/python
from __future__ import print_function
from collections import Counter
""" Prepare histogram of IES cluster sizes by species."""
hist = Counter()
clusterSet = set()
#read silix output line by line
f = open('/home/dsellis/data/IES/analysis/mies/silixout/ies.silixout','r')
for line in f:
    line = line.rstrip()
    (clusterId, iesId) = line.split()
    sp = iesId[0:3]
    iesId = iesId[4:]
    hist[(clusterId, sp)] += 1
    clusterSet.add(clusterId)
#count IESs in each cluster by species
# print('cluster', 'species', 'number', sep = "\t")
# for i in hist:
#     print(i[0], i[1], hist[i], sep = "\t")

spL = ['ppr', 'pbi', 'pte', 'ppe', 'pse', 'poc', 'ptr', 'pso', 'pca']
print("\t".join(["cluster"] + spL + ["total"]), sep = "\t", end = "\n")
for cl in clusterSet:
    L = []
    for sp in spL:
        L.append(hist[(cl, sp)])
    if(sum(L) > 20):
        L = [cl] + L + [sum(L)]
        print("\t".join(str(x) for x in L), end="\n", sep = "\t")


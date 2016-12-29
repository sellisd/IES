#!/usr/bin/python
from __future__ import print_function
from collections import Counter

inputfile = '/Users/dsellis/data/IES/analysis/mies/blastout/mies.blastout'
outputfile = '/Users/dsellis/data/IES/analysis/mies/blastout/mies.blastout.table'

fin = open(inputfile, 'r')
fout = open(outputfile, 'w')
hist = Counter()
clusterSet = set()

for line in fin:
	line = line.rstrip()
	(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
	cluster = qseqid
	species = sseqid[0:3]
	hist[(cluster, species)] += 1
	clusterSet.add(cluster)


spL = ['ppr', 'pbi', 'pte', 'ppe', 'pse', 'poc', 'ptr', 'pso', 'pca']
fout.write("\t".join(["cluster"] + spL + ["total"]) + "\n")
for cl in clusterSet:
    L = []
    for sp in spL:
        L.append(hist[(cl, sp)])
    L = [cl] + L + [sum(L)]
    fout.write("\t".join(str(x) for x in L) + "\n")

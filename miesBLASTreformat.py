#!/usr/bin/python
from __future__ import print_function
from collections import Counter
import sys, getopt

# Prepare histogram of IES cluster sizes by species from BLAST results

inputfile = '';  # blastout/mies.blastout
outputfile = ''; # blastout/mies.blastout.table

usage = "./miesBLASTreformat.py -i <inputfile> -o <outputfile>"

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt == "-i":
        inputfile = arg
    elif opt == "-o":
        outputfile = arg

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

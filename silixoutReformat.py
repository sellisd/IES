#!/usr/bin/python
from __future__ import print_function
from collections import Counter
import sys, getopt
""" Prepare histogram of IES cluster sizes by species."""

inputfile = ''
outputfile = ''
usage = "./silixoutReformat.py -i <inputfile> -o <outputfile> -c cutoff"
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:c:",["ifile=","ofile=","cutoff="])
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-i", "--ifile"):
        inputfile = arg
    elif opt in ("-o", "--ofile"):
        outputfile = arg
    elif opt in ("-c", "--cutoff"):
        cutoff = int(arg)

hist = Counter()
clusterSet = set()

#read silix output line by line
# input /home/dsellis/data/IES/analysis/mies/silixout/ies.silixout

fin = open(inputfile,'r')
fout = open(outputfile, 'w')
for line in fin:
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
fout.write("\t".join(["cluster"] + spL + ["total"]) + "\n")
for cl in clusterSet:
    L = []
    for sp in spL:
        L.append(hist[(cl, sp)])
    if(sum(L) > cutoff):
        L = [cl] + L + [sum(L)]
        fout.write("\t".join(str(x) for x in L) + "\n")


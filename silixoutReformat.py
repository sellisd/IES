#!/usr/bin/python
from __future__ import print_function
from collections import Counter
import sys, getopt

# Prepare histogram of IES cluster sizes by species.

inputfile = ''
outputfile = ''
cutoff = 0
spcutof = 0
usage = "./silixoutReformat.py -i <inputfile> -o <outputfile> -c seqcutoff -s speciesNoCutoff"
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:c:s:",["ifile=","ofile=","cutoff=", "sp="])
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
    elif opt in ("-s", "--sp"):
        spcutof = int(arg)

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
    spNo = 0
    for sp in spL:
        L.append(hist[(cl, sp)])
        if(hist[(cl,sp)] != 0):
            spNo += 1
    if(sum(L) > cutoff and spNo>spcutof):
        L = [cl] + L + [sum(L)]
        fout.write("\t".join(str(x) for x in L) + "\n")


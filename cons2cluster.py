#!/usr/bin/python
from __future__ import print_function
from collections import Counter

# Match BLAST results of consensus against clusters.

silixout = '/home/dsellis/data/IES/analysis/mies/silixout/ies.silixout'
datafile = '/home/dsellis/data/IES/analysis/mies/miescons.dat'
bloutcons = '/home/dsellis/data/IES/analysis/mies/blastout/cons.blastout'

fs = open(silixout, 'r')
fd = open(datafile, 'r')
fb = open(bloutcons, 'r')

iess = {}
# read IES-cluster asignment
for line in fs:
    line = line.rstrip()
    (cluster, ies) = line.split()
    iess[ies] = cluster

# read cluster-consensus correlation
consIds = {}
header = fd.readline()
for line in fd:
    line = line.rstrip()
    if line:
        (cluster, consensusId, peak, lengthL, lengthU) = line.split()
        consIds[consensusId] = cluster

net = Counter()
for line in fb:
    line = line.rstrip()
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
    if iess[sseqid] != consIds[qseqid]:
        net[(consIds[qseqid], iess[sseqid])] += 1

print('from', 'to', 'weights', sep = "\t", end = "\n")
for i in net:
    if(net[i] > 2):
        print(i[0], i[1], net[i], sep = "\t", end = "\n")


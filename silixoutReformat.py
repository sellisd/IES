#!/usr/bin/python
from __future__ import print_function
from collections import Counter
""" Prepare histogram of IES cluster sizes by species."""
hist = Counter()
#read silix output line by line
f = open('/home/dsellis/data/IES/analysis/mies/silixout/ies.silixout','r')
for line in f:
    line = line.rstrip()
    (clusterId, iesId) = line.split()
    sp = iesId[0:3]
    iesId = iesId[4:]
    hist[(clusterId, sp)] += 1
    
#count IESs in each cluster by species
print('cluster', 'species', 'number', sep = "\t")
for i in hist:
    print(i[0], i[1], hist[i], sep = "\t")


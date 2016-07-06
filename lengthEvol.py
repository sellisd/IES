#!/usr/bin/python
from __future__ import print_function
import itertools

"""Length difference in homologous IESs.

    Calculates the absolute difference of IES length (nt) in all pairs of IESs
    within a group of homologous IES."""

alF = '/home/dsellis/data/IES/analysis/tables/ageLength.dat'
elF = '/home/dsellis/data/IES/analysis/tables/lengthEvol.dat'
f = open(alF, 'r')
f.readline() # header
d = {}
for line in f:
    line = line.rstrip()
    (abr, iesId, length, homIESId, geneFamily, age) = line.split("\t")
    d.setdefault((homIESId), []).append(length)

f.close()

o = open(elF, 'w')
for i in d:
    for pair in itertools.combinations(d[i], 2):
        o.write(str(abs(int(pair[0]) - int(pair[1]))) + '\n')
o.close()

#!/usr/bin/python
from __future__ import print_function
from Bio import motifs
from Bio.Seq import Seq
import re
import os

"""
Calculate position weight matrices for IES groups
"""

base = '/home/dsellis/data/IES/analysis/iesdb/'
speciesA = ('ppr', 'pbi', 'pte', 'ppe', 'pse', 'poc', 'ptr', 'pso', 'pca')
baseOut = '/home/dsellis/data/IES/analysis/figures/wlogo'

lengthMin = [18, 34, 41] + [52+n*10 for n in range(0, 44)]
lengthMax = [27, 38, 45] + [56+n*10 for n in range(0, 44)]
lb = range(0, (len(lengthMin)-1))

if not os.path.exists(baseOut):
    os.makedirs(baseOut)

for sp in speciesA:
    fname = os.path.join(base, sp + '.iesdb')
    for lbi in lb: # for each length bin
        f = open(fname, "r")
        f.readline() # header
        fbA = []
        for line in f:
            line.strip()
            (id, scaffold, altSeqNo, startLocs, endLocs, upstreamFlank, downstreamFlank,length, isFloating, merged, inCDS, inInter, inIntron, front, back, seq) = line.split("\t")
            l = int(length)
            if isFloating == '0' and merged == '0' and l >= lengthMin[lbi] and l < lengthMax[lbi]:
                # exclude if have an N
                if re.search('N', front):
                    continue
                if re.search('N', back):
                    continue
                fbA.append(Seq(front + back))
        f.close()

        fbm = motifs.create(fbA)
        peakNo = str("%0.3d" %lbi) 
        figName = os.path.join(baseOut, sp + '.' + peakNo + '.png')
        fbm.weblogo(figName, logo_title = sp + ': ' + peakNo)


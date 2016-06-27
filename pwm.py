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
startX = [18, 34, 41] + [52+n*10 for n in range(0, 20)]
lengthMin = startX[0:-1]
lengthMax = [i-1 for i in startX[1:]]
lengthMin.append(250)
lengthMax.append(5000)
lb = range(0, len(lengthMin))
#lengthMin = [18, 75, 122, 220, 250]
#lengthMax = [33, 81, 131, 240, 10000]
#lb = ['1', '5', '10', '~20', '>250']
if not os.path.exists(baseOut):
    os.makedirs(baseOut)

print( 'species peak min max consensus')

for sp in speciesA:
    fname = os.path.join(base, sp + '.iesdb')
    for lbi in range(0, len(lb)): # for each length bin
        f = open(fname, "r")
        peakNo = str("%0.3d" %lbi)
        figName = os.path.join(baseOut, sp + '.' + peakNo + '.png')
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
                fbA.append(Seq(front[0:5]))
        f.close()

        fbm = motifs.create(fbA)
        title = sp + ': ' + str(lb[lbi]) + '[' + str(lengthMin[lbi]) + '-' + str(lengthMax[lbi]) + ']'
        print( sp + ' ' + str(lb[lbi]) + ' ' + str(lengthMin[lbi]) + ' ' + str(lengthMax[lbi]) + ' ' + fbm.consensus)
        if os.path.exists(figName):
            continue # do not redo figures
#        print(title)
#        fbm.weblogo(figName, logo_title = title)
#        quit()

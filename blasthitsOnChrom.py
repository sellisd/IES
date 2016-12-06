#!/usr/bin/python
from __future__ import print_function
from __future__ import division

# read scaffold lenghts
pprs = open('/home/dsellis/data/IES/analysis/filtscaf/ppr.scaf', 'r')
pbis = open('/home/dsellis/data/IES/analysis/filtscaf/pbi.scaf', 'r')
ptes = open('/home/dsellis/data/IES/analysis/filtscaf/pte.scaf', 'r')
ppes = open('/home/dsellis/data/IES/analysis/filtscaf/ppe.scaf', 'r')
pses = open('/home/dsellis/data/IES/analysis/filtscaf/pse.scaf', 'r')
pocs = open('/home/dsellis/data/IES/analysis/filtscaf/poc.scaf', 'r')
ptrs = open('/home/dsellis/data/IES/analysis/filtscaf/ptr.scaf', 'r')
psos = open('/home/dsellis/data/IES/analysis/filtscaf/pso.scaf', 'r')
pcas = open('/home/dsellis/data/IES/analysis/filtscaf/pca.scaf', 'r')

def readScafLength(f, abr, d):
    f.readline() #header
    for line in f:
        line = line.rstrip()
        (scaffold, length, gc) = line.split()
        d[abr + '.' + scaffold] = int(length)

scafL = {} # scaffold lengths

readScafLength(pprs, "ppr", scafL)
readScafLength(pbis, "pbi", scafL)
readScafLength(ptes, "pte", scafL)
readScafLength(ppes, "ppe", scafL)
readScafLength(pses, "pse", scafL)
readScafLength(pocs, "poc", scafL)
readScafLength(ptrs, "ptr", scafL)
readScafLength(psos, "pso", scafL)
readScafLength(pcas, "pca", scafL)

pprs.close()
pbis.close()
ptes.close()
ppes.close()
pses.close()
pocs.close()
ptrs.close()
psos.close()
pcas.close()

# read blast output
blastf = '/home/dsellis/data/IES/analysis/mies/blastout/miesMac'

f = open(blastf, 'r')
nrow = 1
for line in f:
    line = line.rstrip()
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
    # calculate ratio of hit (start)
    print(str(nrow) + ' ' + str(int(sstart)/scafL[sseqid]))
    nrow += 1

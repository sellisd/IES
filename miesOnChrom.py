#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from pyies.functions import loadScafLengths
import os.path

# Read blast output of mobile IES consensus against all IESs and calculate their location on scaffolds
# default parameters
basePath = '/Users/dsellis/data/IES'
blastf = os.path.join(basePath, "analysis/mies/blastout/mies.blastout")
scafL = {} # scaffold lengths

loadScafLengths(basePath, scafL)

# find in which chromosome IESs are

def miesInScaf(f, abr, d):
    """ Create dictionary with ies as key file with scaffold lengths into dictionary."""
    f.readline() #header
    for line in f:
        line = line.rstrip()
        (iesId,	scaffold, altSeqNo, startLocks) = line.split()[:4]
        iesLoc = int(startLocks.split(",")[0])
        scafLength = scafL[abr + '.' + scaffold]
        d[abr + '.' + iesId] = iesLoc / scafLength

#Load ies in scaffold information from all species.
iesScaf = {}
pprs = open(os.path.join(basePath, 'analysis/iesdb/ppr.iesdb'), 'r')
pbis = open(os.path.join(basePath, 'analysis/iesdb/pbi.iesdb'), 'r')
ptes = open(os.path.join(basePath, 'analysis/iesdb/pte.iesdb'), 'r')
ppes = open(os.path.join(basePath, 'analysis/iesdb/ppe.iesdb'), 'r')
pses = open(os.path.join(basePath, 'analysis/iesdb/pse.iesdb'), 'r')
pocs = open(os.path.join(basePath, 'analysis/iesdb/poc.iesdb'), 'r')
ptrs = open(os.path.join(basePath, 'analysis/iesdb/ptr.iesdb'), 'r')
psos = open(os.path.join(basePath, 'analysis/iesdb/pso.iesdb'), 'r')
pcas = open(os.path.join(basePath, 'analysis/iesdb/pca.iesdb'), 'r')

miesInScaf(pprs, "ppr", iesScaf)
miesInScaf(pbis, "pbi", iesScaf)
miesInScaf(ptes, "pte", iesScaf)
miesInScaf(ppes, "ppe", iesScaf)
miesInScaf(pses, "pse", iesScaf)
miesInScaf(pocs, "poc", iesScaf)
miesInScaf(ptrs, "ptr", iesScaf)
miesInScaf(psos, "pso", iesScaf)
miesInScaf(pcas, "pca", iesScaf)

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
f = open(blastf, 'r')
nrow = 1
for line in f:
    line = line.rstrip()
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
    # calculate ratio of hit (start)
    print(str(nrow) + ' ' + str(iesScaf[sseqid]))
    nrow += 1


#!/usr/bin/python
from __future__ import print_function
from collections import defaultdict

bloutF = '/home/dsellis/data/IES/analysis/mies/blastout/ies.blastout'
flankBlastF = '/home/dsellis/data/IES/analysis/mies/blastout/iesflanks.blastout'
homOutF = '/home/dsellis/data/IES/analysis/mies/blastout/ies.homl.blastout'
nhomOutF = '/home/dsellis/data/IES/analysis/mies/blastout/ies.nonhoml.blastout'

# criteria for homology of flanking sequences
evalueMin = 10^8 # less than evalue
lengthMin = 150  # more or equal to
pidentMin = 75   # more or equal to
# window size for HPS to be considered close to the boundary of an IES
window = 20

# read flank hits
fl = open(flankBlastF, 'r')
hom = {}
for line in fl:
    line = line.rstrip()
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
    if qseqid != sseqid and  float(evalue) < evalueMin and int(length) >= lengthMin and pident >= pidentMin:
        hom[qseqid, sseqid] = 1

fl.close()

#load IES lengths
ppr = open("/home/dsellis/data/IES/analysis/iesdb/ppr.iesdb", 'r')
pbi = open("/home/dsellis/data/IES/analysis/iesdb/pbi.iesdb", 'r')
pte = open("/home/dsellis/data/IES/analysis/iesdb/pte.iesdb", 'r')
ppe = open("/home/dsellis/data/IES/analysis/iesdb/ppe.iesdb", 'r')
pse = open("/home/dsellis/data/IES/analysis/iesdb/pse.iesdb", 'r')
poc = open("/home/dsellis/data/IES/analysis/iesdb/poc.iesdb", 'r')
ptr = open("/home/dsellis/data/IES/analysis/iesdb/ptr.iesdb", 'r')
pso = open("/home/dsellis/data/IES/analysis/iesdb/pso.iesdb", 'r')
pca = open("/home/dsellis/data/IES/analysis/iesdb/pca.iesdb", 'r')

def readIESlength(f, abr, d):
    f.readline() # header
    for line in f:
        line = line.rstrip()
        (iesid, scaffold, altSeqNo, startLocs, endLocs, upstreamFlank, downstreamFlank, length, isFloating, merged, inCDS, inInter, inIntron, front, back, seq) = line.split()
        d[abr + '.' + iesid] = int(length) + 2 # length of sequence including the TAs


iesLengths = {} # IES lengths

readIESlength(ppr, "ppr", iesLengths)
readIESlength(pbi, "pbi", iesLengths)
readIESlength(pte, "pte", iesLengths)
readIESlength(ppe, "ppe", iesLengths)
readIESlength(pse, "pse", iesLengths)
readIESlength(poc, "poc", iesLengths)
readIESlength(ptr, "ptr", iesLengths)
readIESlength(pso, "pso", iesLengths)
readIESlength(pca, "pca", iesLengths)

ppr.close()
pbi.close()
pte.close()
ppe.close()
pse.close()
poc.close()
ptr.close()
pso.close()
pca.close()

hsps = defaultdict(lambda: [False, False, False, False])

bl = open(bloutF, 'r')

for line in bl:
    line = line.rstrip()
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
    qb = int(qstart)
    qe = int(qend)
    sb = int(sstart)
    se = int(send)
    lq = iesLengths[qseqid]
    ls = iesLengths[sseqid]
    if qseqid == sseqid:
        continue # skip self-hits
    if qb <= window and sb <= window:
        hsps[(qseqid, sseqid)][0] = True
    if qe >= lq - window + 1 and se >= ls - window + 1:
        hsps[(qseqid, sseqid)][1] = True
    if qb <= window and sb >= ls - window + 1:
        hsps[(qseqid, sseqid)][2] = True
    if qe >= lq - window + 1 and se <= window:
        hsps[(qseqid, sseqid)][3] = True

bl.seek(0) #rewind

homOut = open(homOutF, 'w')
nhomOut = open(nhomOutF, 'w')

for line in bl:
    line = line.rstrip()
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
    if ((qseqid, sseqid) in hsps):
        print(qseqid+','+sseqid +':' +str(hsps[(qseqid, sseqid)]))
        if hsps[(qseqid, sseqid)][0:2] == [True, True] or hsps[(qseqid, sseqid)][2:] == [True, True]:
            if (qseqid,sseqid) in hom:
                homOut.write(line+"\n")
            else:
                nhomOut.write(line+"\n")

homOut.close()
nhomOut.close()
bl.close()

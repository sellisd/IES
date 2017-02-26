#!/usr/bin/env python
from __future__ import print_function
from collections import Counter
from userOptions import basePath
import os.path, sys, getopt

# Count mobile IES

infile  = os.path.join(basePath, 'analysis', 'mies', 'blastout', 'mies.blastout')
outfile = ''

usage = "./countMobileIES.py -i <inputfile> -o <outputfile>"

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(0)
    elif opt in '-i':
        infile = arg
    elif opt in '-o':
        outfile = arg

mies = Counter()
speciesAbr = ('ppr', 'pbi', 'pte', 'ppe', 'pse', 'poc', 'ptr', 'pso', 'pca')
miesId = set()
with open(infile, 'r') as f:
    for line in f:
        line = line.rstrip()
        (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
        abr = sseqid[0:3]
        mies[(qseqid, abr)] += 1
        miesId.add(qseqid)

with open(outfile, 'w') as fout:
    fout.write('mies' + "\t".join(speciesAbr) + '\n')
    for mi in miesId:
        fout.write(str(mi) + "\t")
        fout.write("\t".join([str(mies[(mi, abr)]) for abr in speciesAbr]))
        fout.write("\n")

#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from pyies.functions import loadScafLengths
import sys, getopt
import os.path

# default parameters
basePath = '/Users/dsellis/data/IES'
#blastf = ""
usage = "./blasthitsOnChrom.py -b <basePath>"

try:
    opts, args = getopt.getopt(sys.argv[1:],"hb:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == "-b":
        basePath = arg

blastf = os.path.join(basePath, "analysis/mies/blastout/miesMac")
scafL = {} # scaffold lengths

loadScafLengths(basePath, scafL)

# read blast output
f = open(blastf, 'r')
nrow = 1
for line in f:
    line = line.rstrip()
    (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = line.split()
    # calculate ratio of hit (start)
    print(str(nrow) + ' ' + str(int(sstart)/scafL[sseqid]))
    nrow += 1

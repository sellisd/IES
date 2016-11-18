#!/usr/bin/python
from __future__ import print_function
from Bio.Seq import Seq
from Bio import SeqIO
import sys, getopt

# default parameters
fl = 100 # flank length
gap = 10
genomeFile = "/home/dsellis/data/IES/genomicData/primaurelia/genome/pprimaurelia_mac_AZ9-3_v1.0.fa"
iesdbFile = "/home/dsellis/data/IES/analysis/iesdb/ppr.iesdb"
outFile = "/home/dsellis/data/IES/analysis/mies/fasta/ppr.fl.fa"
prefix = "ppr" #prepend to IES names

usage = "./flankFasta.py -m <MacGenomeFile> -i <iesdbFile> -o <outputFile> -f flankLength -g gap -p prefix"

try:
    opts, args = getopt.getopt(sys.argv[1:],"hm:i:o:f:g:p:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == "-m":
        genomeFile = arg
    elif opt == "-i":
        iesdbFile = arg
    elif opt == "-o":
        outFile = arg
    elif opt == "-f":
        fl = int(arg)
    elif opt == "-g":
        gap = int(arg)
    elif opt == "-p":
        prefix = arg

# build a dictionary of scaffold lengths and sequences
scafLengths = {}
scafO = {}
for seq_record in SeqIO.parse(genomeFile, "fasta"):
   scafLengths[seq_record.id] = len(seq_record)
   scafO[seq_record.id] = seq_record

# read iesdb
F = [] # list of flanks
f = open(iesdbFile, "r")
f.readline() #header
edgeIESs = 1
for line in f:
    line = line.rstrip()
    (iesId, scaffold, altSeqNo, startLocs, endLocs, upstreamFlank, downstreamFlank, length, isFloating, merged, inCDS, inInter,inIntron, front, back, seq) = line.split("\t")
    b = int(startLocs.split(",")[0]) # keep the first element
    e = int(endLocs.split(",")[0]) # keep the first element
    flank1b = b - gap - fl + 1
    flank1e = b - gap + 1
    flank2b = e + gap
    flank2e = e + gap + fl
    if flank1b < 1 or flank2e > scafLengths[scaffold]:
        edgeIESs += 1
    else:
        flank1 = scafO[scaffold][flank1b - 1 : flank1e - 1]
        flank2 = scafO[scaffold][flank2b - 1 : flank2e - 1]
        merged = flank1 + flank2
        merged.id = prefix + "." + iesId
        merged.description = "merged flanks"
        F.append(merged)

# save output
fout = open(outFile, "w")
SeqIO.write(F, fout, "fasta")
print(prefix + "skipped IESs: " + str(edgeIESs))

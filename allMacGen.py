#!/usr/bin/python
from __future__ import print_function
from Bio.Seq import Seq
from Bio import SeqIO
import os

notationF = '/home/dsellis/data/IES/analysis/notation.csv'
outF = '/home/dsellis/data/IES/tempdat/allMac.fa'
o = open(outF, 'w')
f = open(notationF, 'r')
f.readline() #header
for line in f:
    line = line.rstrip()
    (abbreviation, datapath, binomial, taxId, geneGff, cdsF, protF, geneF, MacF, iesGff, annotation, prefix, micbam) = line.split("\t")[0:13]
    for record in SeqIO.parse(os.path.join(datapath,MacF), "fasta"):
        record.id = abbreviation + '.' + record.id
        SeqIO.write(record, o, "fasta")

# read notation
# for each species read Mac genome
# save all in one file

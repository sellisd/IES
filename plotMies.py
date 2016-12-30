#!/usr/bin/python
from __future__ import print_function
from pyies.functions import speciesAbr2bn
import os.path
from collections import Counter
from ete3 import Tree, TextFace
import sys, getopt

#plot tree with mobile IESs on tips
miesDistrF = "/Users/dsellis/data/IES/analysis/mies/miesDistr.tab"
treeF = "/Users/dsellis/data/IES/analysis/phyldogT2/concatSimple.r.tre"
outgroupName = "Tetrahymena_thermophila"
outputFile = ""

usage = """
Draws a tree next to a table of presence absence data for mobile IES families
Usage:

  ./plotMies.py [OPTIONS]

where OPTIONS can be any of the following:
 -t:    tree file
 -m:    mobile IES presence file
 -o:    outputFile
 -h:    this help screen
"""

try:
    opts, args = getopt.getopt(sys.argv[1:],"ht:m:o:h",)
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt == "-m":
        miesDistrF = arg
    elif opt == "-t":
        treeF = arg
    elif opt == "-o":
        outputFile = arg


t = Tree(treeF)
t.set_outgroup(t&outgroupName)

f = open(miesDistrF, 'r')
f.readline()
mies = {}
miesTotal = Counter()
for line in f:
	line.rstrip()
	(miesFamily, abr, numberOfMembers) = line.split("\t")
	bn = speciesAbr2bn(abr, "full", "_")
	mies[(bn, miesFamily)] = int(numberOfMembers)
	miesTotal[miesFamily] += int(numberOfMembers)

for (column, family) in enumerate(sorted(miesTotal, key = miesTotal.get,reverse = True)):
	for leaf in t:
		if (leaf.name, family) in mies:
			members = mies[(leaf.name, family)]
			if members > 0:
				bgcol = "red"
			else:
				bgcol = "white"
			ts = TextFace(members)
			ts.background.color = bgcol
			leaf.add_face(ts, column = column, position = "aligned")

print(outputFile)
if outputFile:
    t.render(outputFile)
else:
    t.show()
sys.exit(0)

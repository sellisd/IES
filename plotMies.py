#!/usr/bin/python
from __future__ import print_function
from pyies.functions import speciesAbr2bn
import os.path
from collections import Counter
from ete3 import Tree, TextFace, NodeStyle

basePath = "/Users/dsellis/data/IES/"
#plot tree with mobile IESs on tips

miesDistrF = "analysis/mies/miesDistr.tab"
treeF = "analysis/phyldogT2/concatSimple.r.tre"
outgroupName = "Tetrahymena_thermophila"
family = "3214"

t = Tree(os.path.join(basePath, treeF))
t.set_outgroup(t&outgroupName)

f = open(os.path.join(basePath, miesDistrF), 'r')
f.readline()
mies = {}
miesTotal = Counter()
for line in f:
	line.rstrip()
	(miesFamily, abr, numberOfMembers) = line.split("\t")
	bn = speciesAbr2bn(abr, "full", "_")
	mies[(bn, miesFamily)] = numberOfMembers
	miesTotal[miesFamily] += int(numberOfMembers)

for (column, family) in enumerate(sorted(miesTotal, key = miesTotal.get,reverse = True)):
	for leaf in t:
		if (leaf.name, family) in mies:
			leaf.add_face(TextFace(mies[(leaf.name, family)]), column = column, position = "aligned")

t.show()
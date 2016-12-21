#!/usr/bin/python
from __future__ import division
from __future__ import print_function
from ete3 import Tree, TextFace
from pyies.functions import phyldogSpeciesTree
import sys, getopt
import os.path

#/phyldogEventSum.py pathToPhyldog/results/*_Events.txt > ~/data/IES/analysis/phyldogT1/geneDuplLossT2.dat
analysis = '2'
basePath = "/Users/dsellis/data/IES/"
outputFile = ""
scale = 1.0
normalize = 0

usage = """
usage:

./phyldogEventPlot.py [OPTIONS]

    where OPTIONS can be any of the following:
    -b <basePath> Default /Users/dsellis/data/IES
    -a [1|2|3] Choice of species tree analysis to use (Default: 2)
    -o outputFile to render image, if omitted visualize in GUI
	-s float to scale output numbers (default 1)
    -n if provided normalize by branch length, else print numbers
	-h this help screen

The coding for species tree analysis is:
1. PHYLDOG species tree
2. concatenate with simple evolution model
3. concatenate with best codon models

Example:
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -s 0.001 -a 1 -o ~/data/IES/analysis/figures/geneConvLoss1.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -s 0.001 -a 2 -o ~/data/IES/analysis/figures/geneConvLoss2.png &
~/anaconda_ete/bin/python ./phyldogEventPlot.py -n -s 0.001 -a 3 -o ~/data/IES/analysis/figures/geneConvLoss3.png &
"""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hb:a:o:ns:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == "-b":
        basePath = arg
    elif opt == "-a":
        if arg in ['1','2','3']:
            analysis = arg
        else:
            print("Unknown value for option analysis (-a)")
            print(usage)
            sys.exit(1)
    elif opt == "-o":
    	outputFile = arg
    elif opt == "-n":
    	normalize = 1
    elif opt == "-s":
    	scale = float(arg)

brlenTreeFile = ""
if analysis == '1':
	brlenTreeFile  = os.path.join(basePath, "analysis/sgf/topoConstrSimple.treefile")
elif analysis == '2':
	brlenTreeFile = os.path.join(basePath, "analysis/phyldogT2/concatSimple.r.tre")
elif analysis == '3':
	brlenTreeFile = os.path.join(basePath, "analysis/phyldogT3/concat.r.tre")
else:
	sys.exit(1)

eventsfile = os.path.join(basePath, 'analysis/phyldogT1/geneDuplLossT'+ analysis +'.dat') # all geneDuplLossT files are in phyldogT1
phyldogTreeFile = os.path.join(basePath, 'analysis/phyldogT' + analysis + '/results/OutputSpeciesTree_ConsensusNumbered.tree')
outgroupName = "Tetrahymena_thermophila"

t = phyldogSpeciesTree(phyldogTreeFile, brlenTreeFile, outgroupName)

brattr = {}
f = open(eventsfile, 'r')
f.readline() #get rid of header
for line in f:
	line = line.rstrip()
	(speciesNode, duplications, losses) = line.split("\t")
	brattr[speciesNode] = (duplications, losses)

for node in t.traverse():
	if node.PHYLDOGid in brattr:
		dEvents = brattr[node.PHYLDOGid][0] # duplication or gene conversion events
		lEvents = brattr[node.PHYLDOGid][1]
		if node.dist > 10**-10:
			dRate = scale * int(dEvents) / node.dist
			lRate = scale * int(lEvents) / node.dist
			dRate = "+%d" % dRate # keep only integer part
			lRate = "-%d" % lRate
		else:
			dRate = 'NA'
			lRate = 'NA'
		node.add_face(TextFace(dRate), column = 0, position = "branch-top")
		node.add_face(TextFace(lRate), column = 0, position = "branch-bottom")

if outputFile:
    t.render(outputFile)
else:
    t.show()
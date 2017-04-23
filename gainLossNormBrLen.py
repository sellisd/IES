#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import pandas as pd
from ete3 import Tree

#Normalize by branch lengths
spTreeF            = ""
gainLossSumF       = ""
outputF            = ""

usage = """
usage:

./gainLossNormBrLen.py [OPTIONS]

where OPTIONS can be any of the following:

    -g: gainLossSum file
    -t: species tree file with branch lengths
    -o: output file name
    -h: this help screen
""";

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:t:o:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == '-g':
        gainLossSumF = arg
    elif opt == '-t':
        spTreeF = arg
    elif opt == '-o':
        outputF = arg
    else:
        print("Unknown argument" + arg)
        print(usage)
        sys.exit(1)

glsF = pd.read_table(gainLossSumF, sep = "\t", index_col = False)

# normalize by branch length
obspgain = {}
obsploss = {}
t = Tree(spTreeF)
# Sum pgij for each i j and
# sum kgij*ng for each ij

for node in t.iter_descendants():
    fromNode = node.up.ND
    toNode = node.ND
    rowIndexes = (glsF['fromNode']==fromNode) & (glsF['toNode']==toNode)
    obspgain[k] = glsF.loc[rowIndexes,'pcijGain'].sum()
    obsploss[k] = glsF.loc[rowIndexes,'pcijLoss'].sum()

for node in t.iter_descendants():
    fromNode = node.up.ND
    toNode = node.ND
    obspgain[k] /= float(node.dist)
    obsploss[k] /= float(node.dist)

with open(outputF, 'w') as fout:
    fout.write("\t".join(["fromNode", "toNode", "pgain", "ploss\n"]))

    for k in obspgain:
        fout.write("\t".join([str(k[0]), str(k[1]), str(obspgain[k]), str(obsploss[k]) + "\n"]))

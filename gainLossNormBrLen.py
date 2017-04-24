#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import pandas as pd
from ete3 import Tree
from pyies.userOptions import basePath
import os.path
import sys, getopt

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
t = Tree(spTreeF)
# Sum pgij for each i j and
# sum kgij*ng for each ij

# build node pairs and br lenght
branches = {}
for node in t.iter_descendants():
    fromNode = node.up.ND
    toNode = node.ND
    branches[(fromNode, toNode)] = float(node.dist)
sumpgijGain = {}
sumpgijLoss = {}
sumNK       = {}
sumIK       = {}
for k in branches:
    rowIndexes  = (glsF['fromNode']==int(k[0])) & (glsF['toNode']==int(k[1]))
    sumpgijGain[k] = glsF.loc[rowIndexes,'pcijGain'].sum()
    sumpgijLoss[k] = glsF.loc[rowIndexes,'pcijLoss'].sum()
    sumNK[k]       = (glsF.loc[rowIndexes, 'ng'] * glsF.loc[rowIndexes, 'kij']).sum()
    # sumNK: sum (ng * kgij)
    sumIK[k]       = (glsF.loc[rowIndexes, 'Ig'] * glsF.loc[rowIndexes, 'kij']).sum()
    # sumIK: sum (Ig kgij)

RGij = {} # rate of gain
RLij = {} # rate of loss
for k in branches:
    RGij[k] = sumpgijGain[k]/(sumNK[k]*branches[k])
    RLij[k] = sumpgijLoss[k]/(sumIK[k]*branches[k])

with open(outputF, 'w') as fout:
    fout.write("\t".join(["fromNode", "toNode", "rgain", "rloss\n"]))

    for k in branches:
        fout.write("\t".join([str(k[0]), str(k[1]), str(RGij[k]), str(RLij[k]) + "\n"]))

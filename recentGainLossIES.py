#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
from pyies.userOptions import basePath
from pyies.NodeDict import NodeDict
from ete3 import Tree
import sys, getopt
import os.path
import numpy as np
# For each gene tree find terminal(leaf) branches
# foreach X.ReconciledTree find terminal branches

# Read average node probabilities

# search node probabilities for leaf branches and select those with change of probability from x = 0.1 to y = 1
# this means it could have been lost in at least one paralog

# print out analysis number genefamily iesColumn

def prevSpec(nodeO):
    """
    find the closest direct ancestor that corresponds to a speciation node
    Args:
    nodeO ete node object
    """
    while nodeO.up:
        nodeO = nodeO.up
        if(nodeO.Ev == 'S'):
            return nodeO
        else:
            pass

outputPath = ""
geneFamilyId = None
analysis = '2'
cutoff = -0.95

usage = """
usage:

./recentGainLossIES.py [OPTIONS]

    where OPTIONS can be any of the following:
    -g string       Gene family Id, if None provided use all gene families
    -a [1|2|3]      Choice of species tree analysis to use (Default: 2)
    -o <outputPath> Path to create output files
    -c float        Cutoff for change of probability (negative for reduction positive for increase)
    -h              This help screen
"""

try:
    opts, args = getopt.getopt(sys.argv[1:],"ha:o:g:c:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == '-g':
        geneFamilyId = arg
    elif opt == '-c':
        cutoff = float(arg)
    elif opt == '-a':
        analysis = arg
    elif opt == "-o":
        outputPath = arg
        if not os.path.isdir(outputPath):
            print("Warning: Not existing path: " + outputPath)
            print(usage)
            quit(1)

# load node dictionary
nodeDictionary = NodeDict(os.path.join(basePath, "analysis", "tables", "nodeDictionary" + analysis + ".dat"))
# load homologous IES information
homIES = pd.read_csv(os.path.join(basePath, "analysis", "iesdb", "homIESdb" + analysis + ".tab"), sep = "\t")
# load ancestral state probabilities
nodeProbsFile = os.path.join(basePath, "analysis", "tables", "avNodeProb" + analysis + ".dat")
nodeProb = pd.read_csv(nodeProbsFile, sep = "\t", dtype = {'cluster':'str',
                                                           'node': 'str',
                                                           'iesColumn': 'str',
                                                            'presence': np.float64})

print("\t".join(['geneFamily', 'iesColumn', 'node']))

if geneFamilyId == None:
    pass
else:
    phyldogTreeFile = os.path.join(basePath, "analysis", "phyldogT" + analysis, "results", geneFamilyId + ".ReconciledTree")
    print(phyldogTreeFile)
    t = Tree(phyldogTreeFile)
    for leaf in t:
        ancestor = prevSpec(leaf)
        ancestorRB = nodeDictionary.phyldog2rb(geneFamilyId, ancestor.ND)
        offspringRB = nodeDictionary.phyldog2rb(geneFamilyId, leaf.ND)
        # the probability of presence on the ancestor node
        aS = nodeProb.loc[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == ancestorRB)]
        # the probability of presence on the offspring node
        oS = nodeProb.loc[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == offspringRB)]
        aS = aS.reset_index()
        oS = oS.reset_index()
        #which homologous IES columns have a change of probability presence larger than the cuttoff
        lost = oS.presence-aS.presence < cutoff
        if lost.any():
            if not aS.iesColumn.equals(oS.iesColumn):
                print("Error! IES column numbers should match!")
            #homIES['geneFamily' == geneFamilyId]
            #geneFamilyId: oS.cluster[lost]
            #iesColumn: oS.iesColumn[lost]
            #node: leaf.name
            print('\t'.join([oS.cluster[lost].item(), oS.iesColumn[lost].item(), leaf.name]))
# from geneFamilyId, iesColumn and node, find IES
# node can give us gene
#The end result would be a table with columns: genefamily, IESColumn, node, where node is the gene name.

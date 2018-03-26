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


geneFamilyId = None
analysis = '2'
cutoff = 0.95
outputFile = "./recentLoss" + analysis + ".dat"
usage = """
usage:

./recentGainLossIES.py [OPTIONS]

    where OPTIONS can be any of the following:
    -g string       Gene family Id, if None provided use all gene families
    -a [1|2|3]      Choice of species tree analysis to use (Default: 2)
    -o outputFile   Path and file name for output (Default: './recentLoss'+ a + '.dat')
    -c float        Absolute cutoff of probability change:
                    if P_offspring - P_ancestral < - cutoff: loss
                    if P_offspring - P_ancestral > cutoff  : gain
                    else                                   : undetermined
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
        outputFile = arg
        if os.path.isfile(outputFile):
            print("Warning: overwriting file " + outputFile + "!")

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


with open(outputFile, 'w') as f:
    f.write("\t".join(['geneFamily', 'iesColumn', 'node\n']))

def printRecentIESLoss(geneFamilyId = 10000, analysis = 2, cutoff = 0.95):
    """Extract the recent IES loss
    Args:
    geneFamilyId:  string  Unique id of gene family
    analysis:      string  Choice of species tree analysis to use (Default: 2)
    -c float        Absolute cutoff of probability change:
                    if P_offspring - P_ancestral < - cutoff: loss
                    if P_offspring - P_ancestral > cutoff  : gain
                    else                                   : undetermined
    """
    phyldogTreeFile = os.path.join(basePath, "analysis", "phyldogT" + analysis, "results", geneFamilyId + ".ReconciledTree")
    t = Tree(phyldogTreeFile)
    for leaf in t:
        ancestor = prevSpec(leaf)
        if ancestor is None:
            # skip leafs that are not preceded by a speciation event
            continue
        ancestorRB = nodeDictionary.phyldog2rb(geneFamilyId, ancestor.ND)
        offspringRB = nodeDictionary.phyldog2rb(geneFamilyId, leaf.ND)
        if (ancestorRB is None) or (offspringRB is None):
            continue
        # the probability of presence on the ancestor node
        aS = nodeProb.loc[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == ancestorRB)]
        # the probability of presence on the offspring node
        oS = nodeProb.loc[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == offspringRB)]
        aS = aS.reset_index()
        oS = oS.reset_index()
        #which homologous IES columns have a change of probability presence larger than the cuttoff
        probabilityChange = oS.presence-aS.presence
        losses = probabilityChange < - cutoff
        gains = probabilityChange > cutoff
        if losses.any() or gains.any():
            if not aS.iesColumn.equals(oS.iesColumn):
                print("Error! IES column numbers should match!")
            for i, dp in enumerate(probabilityChange.tolist()):
                event = None
                if losses[i].item() == True:
                    event = "loss"
                elif gains[i].item() == True:
                    event = "gain"
                else:
                    event = "undetermined"
                with open(outputFile, 'a') as f:
                    f.write('\t'.join([oS.cluster[i], oS.iesColumn[i], leaf.name, event + "\n"]))
            #homIES['geneFamily' == geneFamilyId]
            #geneFamilyId: oS.cluster[lost]
            #iesColumn: oS.iesColumn[lost]
            #node: leaf.name

if geneFamilyId == None:
    #Run the analysis for all gene families for which we have ancestral state reconstructions
    for geneFamilyId in nodeProb.cluster.unique():
        printRecentIESLoss(geneFamilyId, analysis, cutoff)
else:
    printRecentIESLoss(geneFamilyId, analysis, cutoff)

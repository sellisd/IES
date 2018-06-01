#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
from pyies.userOptions import basePath
from pyies.NodeDict import NodeDict
from ete3 import Tree
import sys, getopt
import os.path
import numpy as np

#~/anaconda_ete/bin/python branchEvents.py -n -a 3 -f 6 -t 9 -o branchEvents3_PocPte.dat
#~/anaconda_ete/bin/python branchEvents.py -n -a 3 -f 14 -t 15 -o branchEvents3_PprPpe.dat
# For each gene tree find terminal(leaf) branches
# foreach X.ReconciledTree find terminal branches

# Read average node probabilities

# search node probabilities for leaf branches and select those with change of probability from x = 0.1 to y = 1
# this means it could have been lost in at least one paralog

# print out analysis number genefamily iesColumn

def prevSpecRecursive(nodeO):
    """
    Recursively go up the tree to find the closest direct ancestor that corresponds to a speciation node
    Args:
    nodeO ete node object
    """
    while nodeO.up:
        nodeO = nodeO.up
        if(nodeO.Ev == 'S'):
            return nodeO
        else:
            pass

def prevSpeciation(nodeO):
    """
    Return parent node if it is a speciation node else return None
    Args:
    nodeO ete node object
    """
    ancestor = nodeO.up
    if ancestor:
        if ancestor.Ev == 'S':
            return ancestor
    return None

leafNames = True
recent = False
recursive = False
geneFamilyId = None
analysis = '3'
outputFile = "./branchEvents" + analysis + ".dat"
usage = """
usage:
Classifies and prints events on branches. If option -c is set only the most recent (terminal) branches are considered.
./branchEvents.py [OPTIONS]

    where OPTIONS can be any of the following:
    -g string       Gene family Id, if None provided use all gene families
    -r              Branches can include duplication nodes (recursive search upstream)
                    (By default this is false). This option is not compatible
                    with -f and -t
    -a [1|2|3]      Choice of species tree analysis to use (Default: 3)
    -o outputFile   Path and file name for output (Default: './branchEvents'+ a + '.dat')
    -c              If true only events on the most recent (terminal) branches are considered
    -f string       Upstream species node notation Id for branch definition
    -t string       Downstream species node notation ID for branch definition
    -n              Include an extra column in the output with the name of the leafs
    -h              This help screen
"""

try:
    opts, args = getopt.getopt(sys.argv[1:],"ha:o:g:c:rf:t:n")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == '-r':
        recursive = True
    elif opt == '-g':
        geneFamilyId = arg
    elif opt == '-a':
        analysis = arg
    elif opt == '-c':
        recent = True
    elif opt == '-f':
        fromNode = arg
    elif opt == '-t':
        toNode = arg
    elif opt == '-n':
        leafNames = True
    elif opt == '-o':
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
    outputList = ['geneFamily', 'iesColumn', 'node']
    if leafNames:
        outputList.append("genes")
    outputList[-1] = outputList[-1] + '\n'
    f.write("\t".join(outputList))


def printRecentIESLoss(geneFamilyId = 10000, analysis = 2):
    """Extract the recent IES loss
    Args:
    geneFamilyId:  string  Unique id of gene family
    analysis:      string  Choice of species tree analysis to use (Default: 2)
    """
    phyldogTreeFile = os.path.join(basePath, "analysis", "phyldogT" + analysis, "results", geneFamilyId + ".ReconciledTree")
    t = Tree(phyldogTreeFile)
    for leaf in t.traverse("preorder"):
        ancestor = None
        geneString = ""
        if recursive:
            ancestor = prevSpecRecursive(leaf)
        else:
            ancestor = prevSpeciation(leaf)
        if leafNames:
            #get string with all
            geneNodes = leaf.get_leaves()
            nameList = [node.name for node in geneNodes]
            geneString = ",".join(nameList)
        if ancestor is None:
            # skip leafs that are not preceded by a speciation event
            continue
        ancestorRB = nodeDictionary.phyldog2rb(geneFamilyId, ancestor.ND)
        offspringRB = nodeDictionary.phyldog2rb(geneFamilyId, leaf.ND)
        if (ancestorRB is None) or (offspringRB is None):
            continue
        # apply fiering conditions
        if recent:
            if not leaf.is_leaf():
                continue
        else:
            if ancestor.S != fromNode or leaf.S != toNode:
                continue
        # the probability of presence on the ancestor node
        aS = nodeProb.loc[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == ancestorRB)]
        # the probability of presence on the offspring node
        oS = nodeProb.loc[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == offspringRB)]
        aS = aS.reset_index()
        oS = oS.reset_index()
        #which homologous IES columns have a change of probability presence larger than the cuttoff
        probabilityChange = oS.presence-aS.presence
        classIES = None
        i = 0
        for pAncestral, pOffspring in zip(aS.presence.tolist(), oS.presence.tolist()):
            if pOffspring <= 0.01:
                if pAncestral > 0.95:
                    classIES = "lost"
                elif pAncestral < 0.05:
                    classIES = "absent"
                elif pAncestral > 0.5:
                    classIES = "probably_lost"
                else:
                    classIES = "probably_absent"
            elif pOffspring >= 0.99:
                if pAncestral > 0.95:
                    classIES = "conserved"
                elif pAncestral < 0.05:
                    classIES = "gained"
                elif pAncestral < 0.50:
                    classIES = "probably_gained"
                else:
                    classIES = "probably_conserved"
            else:
                classIES = "undetermined"
            with open(outputFile, 'a') as f:
                f.write('\t'.join([oS.cluster[i], oS.iesColumn[i], leaf.name, classIES, geneString + "\n"]))
            i += 1

if geneFamilyId == None:
    #Run the analysis for all gene families for which we have ancestral state reconstructions
    for geneFamilyId in nodeProb.cluster.unique():
        printRecentIESLoss(geneFamilyId, analysis)
else:
    printRecentIESLoss(geneFamilyId, analysis)

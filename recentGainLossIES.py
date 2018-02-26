#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
from pyies.userOptions import basePath
from pyies.NodeDict import NodeDict
from ete3 import Tree
import sys, getopt
import os.path

# For each gene tree find terminal(leaf) branches
# foreach X.ReconciledTree find terminal branches

# Read average node probabilities

# check if necessary to use the node nodeDictionary

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
cutoff = 1

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
# load ancestral state probabilities
nodeProbsFile = os.path.join(basePath, "analysis", "tables", "avNodeProb" + analysis + ".dat")
nodeProb = pd.read_csv(nodeProbsFile, sep = "\t")

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
        print(ancestorRB + "-" + offspringRB)
        aS = nodeProb.presence[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == ancestorRB)]
        oS = nodeProb.presence[(nodeProb.cluster == geneFamilyId) & (nodeProb.node == offspringRB)]
        oS-aS



print(np.head())
quit()
# output files
rgl = open(os.path.join(outputPath, 'recentGainLossIES.dat'), 'w')

# read node dictionary
rb2phyldog = {}
for line in open(os.path.join(basePath, 'IES/rdb/nodeDictionary.dat')):
    line = line.rstrip()
    (cluster, r, phyldog, rb) = line.split()
    rb2phyldog[(cluster, rb)] = phyldog

# read nodeAsr and make dictionary with cluster node and presence
f = open('/home/dsellis/data/IES_data/msas/nodeAsr.dat')
d={} # dictionary with (cluster, node) = { iesCol1: presence, iesCol2:presence...}
header = f.readline()
for line in f:
    line = line.rstrip()
    (cluster, node, iesColumn, presence) = line.split()
    if(cluster, node) in rb2phyldog:
        nodeP = rb2phyldog[(cluster, node)]
    else:
#        print("skipping :" + cluster + ', ' + node, file = sys.stderr)
        continue
    if (cluster, nodeP) in d:
        d[(cluster, nodeP)][iesColumn] = presence
    else:
        d[(cluster, nodeP)] = {iesColumn : presence}

# for each gene family find terminal branch and walk backwards until the most recent speciation event.
phyldogPath = '/home/dsellis/data/IES_data/msas/phyldog/results/'
inputF = open('/home/dsellis/data/IES_data/msas/asr/geneFamilies.dat', 'r')

clusters = inputF.readlines()
clusters = [i.rstrip() for i in clusters]
headerOut = '\t'.join(['cluster', 'iesColumn', 'gene', 'fromS', 'toS', 'fromNodeP', 'toNodeP', 'presenceFrom', 'presenceTo', 'eventType']) + '\n'
fnt.write(headerOut)
fgl.write(headerOut)
for cluster in clusters:
    fileNameString = phyldogPath+str(cluster)+'.ReconciledTree'
    print(fileNameString)
    t = Tree(fileNameString)
#    print(t.get_ascii(attributes=["Ev","S","ND","name"]))
    for leaf in t:
        anc = prevSpec(leaf)
        if(anc):
            if(cluster, leaf.ND) in d:
                mykeys = [(cluster, anc.ND), (cluster, leaf.ND)]
                probs = [d[x] for x in mykeys]
                # for each ies column
                iesColumns = list(probs[0].keys())
                for i in iesColumns:
                    wasPresent = float(d[(cluster, anc.ND)][i])
                    isPresent = float(d[(cluster, leaf.ND)][i])
                    if isPresent == 0 and wasPresent > 0.99:
                        fgl.write("\t".join([cluster, i, leaf.name, anc.S, leaf.S, anc.ND, leaf.ND, str(wasPresent), str(isPresent), 'loss']) + '\n')
                    elif wasPresent < 0.01 and isPresent == 1:
                        fgl.write("\t".join([cluster, i, leaf.name, anc.S, leaf.S, anc.ND, leaf.ND, str(wasPresent), str(isPresent), 'gain']) + '\n')
                    elif wasPresent < 0.01 and isPresent <0.01:
                        fnt.write("\t".join([cluster, i, leaf.name, anc.S, leaf.S, anc.ND, leaf.ND, str(wasPresent), str(isPresent), 'absent']) + '\n')
                    elif wasPresent > 0.99 and isPresent <0.99:
                        fnt.write("\t".join([cluster, i, leaf.name, anc.S, leaf.S, anc.ND, leaf.ND, str(wasPresent), str(isPresent), 'present']) + '\n')
                    else:
                        fnt.write("\t".join([cluster, i, leaf.name, anc.S, leaf.S, anc.ND, leaf.ND, str(wasPresent), str(isPresent), 'notSign']) + '\n')
        else:
            pass
fgl.close()
fnt.close()
# try:
        #     print(leaf.ND + "\t" + anc.ND)
        # except:
        #     print(cluster + leaf.name)
        #     print(leaf)
        #     print(t.get_ascii(attributes=["Ev","S","ND","name"]))
        #     quit()


# Compare the probability of presence in the two nodes.

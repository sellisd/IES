#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree
import sys

"""
find IES that are lost in a terminal branch
"""
def prevSpec(nodeO):
    """
    find the closest direct ancestor that corresponds to a speciation node
    """
    while nodeO.up:
        nodeO = nodeO.up
        if(nodeO.Ev == 'S'):
            return nodeO
        else:
            pass

# read node dictionary
rb2phyldog = {}
for line in open('/home/dsellis/data/IES_data/rdb/nodeDictionary.dat'):
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
        next
    if (cluster, nodeP) in d:
        d[(cluster, nodeP)][iesColumn] = presence
    else:
        d[(cluster, nodeP)] = {iesColumn : presence}

# for each gene family find terminal branch and walk backwards until the most recent speciation event.
phyldogPath = '/home/dsellis/data/IES_data/msas/phyldog/results/'
inputF = open('/home/dsellis/data/IES_data/msas/asr/geneFamilies.dat', 'r')

clusters = inputF.readlines()
clusters = [i.rstrip() for i in clusters]
print('\t'.join(['cluster', 'iesColumn', 'gene', 'fromS', 'toS', 'fromNodeP', 'toNodeP', 'presenceFrom', 'presenceTo']))
for cluster in clusters:
    fileNameString = phyldogPath+str(cluster)+'.ReconciledTree'
    # print(fileNameString)
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
                        print("\t".join([cluster, i, leaf.name, anc.S, leaf.S, anc.ND, leaf.ND, str(wasPresent), str(isPresent)]))
        else:
            pass

# try:
        #     print(leaf.ND + "\t" + anc.ND)
        # except:
        #     print(cluster + leaf.name)
        #     print(leaf)
        #     print(t.get_ascii(attributes=["Ev","S","ND","name"]))
        #     quit()


# Compare the probability of presence in the two nodes.



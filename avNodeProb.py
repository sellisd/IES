#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from pyies.functions import phyldogSpeciesTree
from pyies.userOptions import basePath
from collections import Counter
from ete3 import Tree, NodeStyle, TextFace
import os.path

# read spEvents into dictionary node->sp event
spEventsF = os.path.join(basePath, 'analysis', 'tables', 'spEvents1.dat')
avNodeProbF = os.path.join(basePath, 'analysis', 'tables', 'avNodeProb1.dat')
nodeDictF = os.path.join(basePath, 'analysis', 'tables', 'nodeDictionary1.dat')
brlenFile = os.path.join(basePath, 'analysis', 'sgf', 'topoConstrSimple.treefile')
phyldogTreeFile = os.path.join(basePath, 'analysis', 'phyldogT1', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree')
outgroupName = "Tetrahymena_thermophila"
outputFile = os.path.join(basePath, 'analysis', 'figures', 'avnodeProb.png')
includeGF = "/Users/dsellis/Desktop/temp.dat"

# if parameter defined load a list of gene families to include from the analysis
includedGeneFamilies = []
if includeGF:
    with open(includeGF, 'r') as f:
        for line in f:
            includedGeneFamilies = [line.rstrip() for line in f]


# load node dictionary
rb2phyldog = {}
with open(nodeDictF, 'r') as f:
    f.readline()
    for line in f:
        line = line.rstrip()
        (cluster, r, phyldog, rb) = line.split("\t")
        rb2phyldog[(cluster, rb)] = phyldog

# load node id speciation events (ND:S)
node2Event = {}
with open(spEventsF, 'r') as f:
    f.readline()
    for line in f:
        line = line.rstrip()
        (cluster, nodeP, spEvent) = line.split("\t")
        node2Event[(cluster, nodeP)] = spEvent

# read avnodProb use dictionary to translate from node to sp event
# sum for each sp event
sumProb = Counter()
countProb = Counter()
with open(avNodeProbF, 'r') as f:
    f.readline()
    for line in f:
        line = line.rstrip()
        (cluster, rb, iesColumn, presence) = line.split("\t")
        #homIESid = column2homIESid[(cluster, iesColumn)]
        nodeP = rb2phyldog[(cluster, rb)]
        if cluster in includedGeneFamilies:
            if (cluster, nodeP) in node2Event: # if node is speciation node (not duplication)
                spEvent = node2Event[(cluster,nodeP)]
                sumProb[spEvent] += float(presence)
                countProb[spEvent] += 1

# load species tree with branch lengths
# load phyldog tree
t = phyldogSpeciesTree(phyldogTreeFile, brlenFile, outgroupName)
for node in t.traverse():
    #print(node.PHYLDOGid+' '+ str(sumProb[node.PHYLDOGid]/countProb[node.PHYLDOGid]))
    #nstyle = NodeStyle()
    p = 100 * sumProb[node.PHYLDOGid]/countProb[node.PHYLDOGid]
    #node.set_style(nstyle)
    node.add_face(TextFace(str(round(p,2))), column = 0, position = "branch-right")

t.render(outputFile)
#for k in countProb:
#    print("\t".join([k, str(sumProb[k]/countProb[k])]))

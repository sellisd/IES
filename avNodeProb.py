#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from pyies.functions import phyldogSpeciesTree
from pyies.userOptions import basePath
from collections import Counter
from ete3 import Tree, NodeStyle, TextFace
import os.path

# read spEvents into dictionary node->sp event
#asrRun = 1|2|3
asrRun = 3
spEventsF = os.path.join(basePath, 'analysis', 'tables', 'spEvents'+str(asrRun)+'.dat')
avNodeProbF = os.path.join(basePath, 'analysis', 'tables', 'avNodeProb'+str(asrRun)+'.dat')
nodeDictF = os.path.join(basePath, 'analysis', 'tables', 'nodeDictionary'+str(asrRun)+'.dat')
if asrRun == 1:
    brlenFile = os.path.join(basePath, 'analysis', 'sgf', 'topoConstrSimple.treefile')
elif asrRun == 2:
    brlenFile = os.path.join(basePath, 'analysis', 'sgf', 'concatSimple.nexus.treefile')
elif asrRun == 3:
    brlenFile = os.path.join(basePath, 'analysis', 'sgf', 'concat.nexus.treefile')
else:
    quit(1)
phyldogTreeFile = os.path.join(basePath, 'analysis', 'phyldogT'+str(asrRun), 'results', 'OutputSpeciesTree_ConsensusNumbered.tree')
outgroupName = "Tetrahymena_thermophila"
outputFile = os.path.join(basePath, 'analysis', 'figures', 'avnodeProb'+str(asrRun)+'.png')
includeGF = "/Volumes/WDC/data/IES/analysis/tables/singleGeneFamilies.dat"

# if parameter defined load a list of gene families to include from the analysis
includedGeneFamilies = []
if includeGF:
    with open(includeGF, 'r') as f:
        includedGeneFamilies = [line.rstrip() for line in f]

#includedGeneFamilies = ['10623']
#print(includedGeneFamilies)

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
        if cluster in includedGeneFamilies:
            #print(line)
            node2Event[(cluster, nodeP)] = spEvent

#quit()
#for k in node2Event:
#    print((k,node2Event[k]))
#quit()
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
                #print(nodeP + ' ' + spEvent+' '+presence)
                sumProb[spEvent] += float(presence)
                countProb[spEvent] += 1.
#quit(0)
# load species tree with branch lengths
# load phyldog tree
t = phyldogSpeciesTree(phyldogTreeFile, brlenFile, outgroupName)
for node in t.traverse():
    #print(node.PHYLDOGid+' '+ str(sumProb[node.PHYLDOGid]/countProb[node.PHYLDOGid]))
    #nstyle = NodeStyle()
    if(countProb[node.PHYLDOGid] == 0):
        p = -1
    else:
        p = 100 * sumProb[node.PHYLDOGid]/countProb[node.PHYLDOGid]
    #node.set_style(nstyle)
    node.add_face(TextFace(str(round(p,2))), column = 0, position = "branch-right")

t.render(outputFile)
#for k in countProb:
#    print("\t".join([k, str(sumProb[k]/countProb[k])]))

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
for asrRun in [1, 2, 3]:
    print(asrRun)
    spEventsF = os.path.join(basePath, 'analysis', 'tables', 'spEvents'+str(asrRun)+'.dat')
    avNodeProbF = os.path.join(basePath, 'analysis', 'tables', 'avNodeProb'+str(asrRun)+'.dat')
    nodeDictF = os.path.join(basePath, 'analysis', 'tables', 'nodeDictionary'+str(asrRun)+'.dat')
    if asrRun == 3:
        brlenFile = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree3b.nhx')
    else:
        brlenFile = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree' + str(asrRun) + '.nhx')
    outputFile = os.path.join(basePath, 'analysis', 'figures', 'avnodeProb'+str(asrRun)+'.png')

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
            nodeP = rb2phyldog[(cluster, rb)]
            if (cluster, nodeP) in node2Event: # if node is speciation node (not duplication)
                spEvent = node2Event[(cluster,nodeP)]
                sumProb[spEvent] += float(presence)
                countProb[spEvent] += 1.

    # load species tree with branch lengths
    t = Tree(brlenFile)
    for node in t.traverse():
        p = 100 * sumProb[node.ND]/countProb[node.ND]
        node.add_face(TextFace(str(round(p,2))), column = 0, position = "branch-right")

    t.render(outputFile)

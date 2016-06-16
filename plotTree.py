#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree, NodeStyle, SeqMotifFace
import pprint
import os.path
import sys
"""
Draw gene family trees with nodes colored by event-type and motif of IESs
"""

pp = pprint.PrettyPrinter(indent=8)

cp = ["#c24c6e", "#ba4c46", "#bf5a2f", "#d39036", "#9b7532", "#b4ad55", "#b6b638", "#88af46", "#477b2f", "#63c471", "#48c19a", "#628ed6", "#6a70d7", "#533585", "#bb85d6", "#b759b8", "#d56cad", "#892863"]
# generated from http://tools.medialab.sciences-po.fr/iwanthue/

class Vividict(dict):
    # by http://stackoverflow.com/users/541136/aaron-hall
    # from http://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries-in-python/19829714#19829714
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

def colorNodes( t, cp ):
    """ Add styling to nodes of a PHYLDOG tree

    t: tree
    cp: list of colors (hex) corresponding to speciation events
    """

    for node in t.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 10
        if node.Ev == 'S':
            nstyle["fgcolor"] = cp[int(node.S) - 1]
            nstyle["shape"] = "circle"
        elif node.Ev == 'D':
            nstyle["fgcolor"] = "black" #cblue4
            nstyle["shape"] = "square"
        else:
            quit()
        node.set_style(nstyle)
    return t


# plot a PHYLDOG tree with annotations for speciation, duplication and alignment/ IES presence/absence
#different colors for each speciation event? or different colors for speciation/duplication

# plot the species tree with annotations for rates of insertion and loss
# sptreeF = '/home/dsellis/data/IES/analysis/phyldog/results/OutputSpeciesTree_ConsensusDuplications.tree'
# st = Tree(sptreeF)
# st.show()
# colorNodes(st, cp)
# st.show()
# quit()

#geneFamily = '1000'
geneFamily = sys.argv[1]

#  load gene family tree
treeF = '/home/dsellis/data/IES/analysis/phyldog/results/' + geneFamily + '.ReconciledTree'
t = Tree(treeF)
if not os.path.isfile(treeF):
    quit("charMatsF is not a file")

# # set styles
t = colorNodes(t, cp)

# load character matrices
charMatsF = '/home/dsellis/data/IES/analysis/iesdb/charMats.tab'
f = open (charMatsF, 'r')
header = f.readline()
charMat = Vividict()
for line in f:
    line = line.rstrip()
    (cluster, column, geneId, begin, end, ies, iesId, beginMSA, endMSA) = line.split("\t")
    charMat[cluster][column][geneId] = [begin, end, ies, iesId, beginMSA, endMSA]
#    charMat[cluster] = {column: {geneId: [begin, end, ies, iesId, beginMSA, endMSA]}}
    
# add motif faces
for l in t:
#    print(l.name)
    # build motif
    iesmotif = []
    for column in charMat[geneFamily]:
#        print(column)
        for geneId in charMat[geneFamily][column]:
            if geneId == l.name:
                begin    = charMat[geneFamily][column][geneId][0]
                end      = charMat[geneFamily][column][geneId][1]
                ies      = charMat[geneFamily][column][geneId][2]
                iesId    = charMat[geneFamily][column][geneId][3]
                beginMSA = charMat[geneFamily][column][geneId][4]
                endMSA   = charMat[geneFamily][column][geneId][5]
                iesmotif.append( [int(begin), int(end), "line", None, 10, "black", "blue", None]) # add range of homologous IES as line
                if iesId != 'NA': # if present add IES as rectangles
                    iesmotif.append( [int(beginMSA), int(endMSA), "[]", None, 10, "black", "red", None] )
    if iesmotif:
        seqFace = SeqMotifFace(seq = None, motifs = iesmotif, gap_format = "line")
        l.add_face(seqFace, 0, "aligned")
#           print("  " + geneId)
# # load alignment add biopython
# #alnF = '/home/dsellis/data/analysis/msas/filtered/cluster.10000.nucl.fa'

t.show()

 

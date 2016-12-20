#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree, SeqMotifFace, ImgFace, TreeStyle, SeqGroup
import os.path
import sys, getopt
import string
from collections import defaultdict
from pyies.functions import *

# plot gene family phylogeny with diagram of IESs in MSA

# parameters
basePath = "/home/dsellis/data/IES/"
geneFamily = ""
phyldogResultsPath = "analysis/phyldogT1/results/"
charMatFile = "analysis/iesdb/charMats.tab"
outputFile = ""

usage = "./plotTree.py -b <basePath> -g <geneFamily> -p <phyldogPath> -m <charMatFile> -o <outputFile>"

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:p:m:o:b:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == "-g":
        geneFamily = arg
    elif opt == "-p":
        phyldogPath = arg
    elif opt == "-m":
        charMatFile = arg
    elif opt == "-o":
        outputFile = arg
    elif opt == "-b":
        basePath = arg

#  load gene family tree
geneFamFile = geneFamily + ".ReconciledTree"
treeF = os.path.join(basePath, phyldogResultsPath, geneFamFile)

t = Tree(treeF)

t = colorNodes(t, 1)
# load character matrices
f = open (os.path.join(basePath, charMatFile), 'r')
header = f.readline()
#charMat = Vividict()
gfGenes  = defaultdict(list) # dictionary key is gene family value list of genes
gfhomIES = defaultdict(list) # dictionary key is gene family value list of homologous IESs
charMat = {}
for line in f:
    line = line.rstrip()
    (cluster, column, geneId, begin, end, ies, iesId, beginMSA, endMSA) = line.split("\t")
    gfGenes[cluster].append(geneId)
    gfhomIES[cluster].append(column)
    charMat[(cluster, column, geneId)] = [begin, end, ies, iesId, beginMSA, endMSA]
#    charMat[cluster][column][geneId] = [begin, end, ies, iesId, beginMSA, endMSA]
#    charMat[cluster] = {column: {geneId: [begin, end, ies, iesId, beginMSA, endMSA]}}

# load mean ancestral states
# asrF = '/home/dsellis/data/IES/analysis/tables/avNodeProb.dat';
# f = open (asrF, 'r')
# header = f.readline()

# load nucleotide sequences for all genes!
nuclAlnFile = os.path.join(basePath, "analysis/msas/filtered/cluster." + geneFamily + ".nucl.fa")
seqs = SeqGroup(sequences = nuclAlnFile, format = "fasta")

for leaf in t:
    geneId = leaf.name
    seq = seqs.get_seq(geneId)
    seq = seq.translate(None, string.ascii_lowercase) # keep only CDS 
    iesmotif = [[1, len(seq), "line", 2, 5, None, None, None]]
    for homIES in gfhomIES[geneFamily]:
        (begin, end, ies, iesId, beginMSA, endMSA) = charMat[(geneFamily, homIES, geneId)]
        if ies == '?':
            iesmotif.append([int(beginMSA), int(endMSA),"()", 10, 10, "red", "black", "arial|8|black|?"])
        elif ies == '1':
            iesmotif.append([int(beginMSA), int(endMSA),"[]", 10, 10, "black", "red", "arial|8|black|" + iesId])
        elif ies == '0':
            iesmotif.append([int(begin), int(end), "[]", 10, 10, "silver", "silver", None]) 
        else:
            quit(1)
    seqFace = SeqMotifFace(seq = seq, motifs = iesmotif, gap_format = "blank", seq_format = "line")
    leaf.add_face(seqFace, 0, "aligned")


#           print("  " + geneId)
# # load alignment add biopython
# #alnF = '/home/dsellis/data/analysis/msas/filtered/cluster.10000.nucl.fa'
if outputFile:
    t.render(outputFile)
else:
    t.show()

"""
t.show()

# Draw trees.

pp = pprint.PrettyPrinter(indent=8)

# SPECIATION TREE #
###################
wgd1 = Tree('((P_caudatum:1[&&NHX:Ev=S:S=3:ND=3],(((P_sexaurelia:1[&&NHX:Ev=S:S=7:ND=7],P_sonneborni:1[&&NHX:Ev=S:S=8:ND=8]):1[&&NHX:Ev=S:S=5:ND=5],(((P_pentaurelia:1[&&NHX:Ev=S:S=13:ND=13],P_primaurelia:1[&&NHX:Ev=S:S=14:ND=14]):1[&&NHX:Ev=S:S=11:ND=11],(P_biaurelia:1[&&NHX:Ev=S:S=15:ND=15],(P_octaurelia:1[&&NHX:Ev=S:S=17:ND=17],P_tetraurelia:1[&&NHX:Ev=S:S=18:ND=18]):1[&&NHX:Ev=S:S=16:ND=16]):1[&&NHX:Ev=S:S=12:ND=12]):1[&&NHX:Ev=S:S=9:ND=9],P_tredecaurelia:1[&&NHX:Ev=S:S=10:ND=10]):1[&&NHX:Ev=S:S=6:ND=6]):1[&&NHX:Ev=S:S=4:ND=4],((P_sexaurelia:1[&&NHX:Ev=S:S=7:ND=7],P_sonneborni:1[&&NHX:Ev=S:S=8:ND=8]):1[&&NHX:Ev=S:S=5:ND=5],(((P_pentaurelia:1[&&NHX:Ev=S:S=13:ND=13],P_primaurelia:1[&&NHX:Ev=S:S=14:ND=14]):1[&&NHX:Ev=S:S=11:ND=11],(P_biaurelia:1[&&NHX:Ev=S:S=15:ND=15],(P_octaurelia:1[&&NHX:Ev=S:S=17:ND=17],P_tetraurelia:1[&&NHX:Ev=S:S=18:ND=18]):1[&&NHX:Ev=S:S=16:ND=16]):1[&&NHX:Ev=S:S=12:ND=12]):1[&&NHX:Ev=S:S=9:ND=9],P_tredecaurelia:1[&&NHX:Ev=S:S=10:ND=10]):1[&&NHX:Ev=S:S=6:ND=6]):1[&&NHX:Ev=S:S=4:ND=4]):1[&&NHX:Ev=D:S=4:ND=4]):1[&&NHX:Ev=S:S=1:ND=1],T_thermophila:1[&&NHX:Ev=S:S=2:ND=2])[&&NHX:Ev=S:S=0:ND=0];')

basest = Tree('((P_caudatum:1[&&NHX:Ev=S:S=3:ND=3],((P_sexaurelia:1[&&NHX:Ev=S:S=7:ND=7],P_sonneborni:1[&&NHX:Ev=S:S=8:ND=8]):1[&&NHX:Ev=S:S=5:ND=5],(((P_pentaurelia:1[&&NHX:Ev=S:S=13:ND=13],P_primaurelia:1[&&NHX:Ev=S:S=14:ND=14]):1[&&NHX:Ev=S:S=11:ND=11],(P_biaurelia:1[&&NHX:Ev=S:S=15:ND=15],(P_octaurelia:1[&&NHX:Ev=S:S=17:ND=17],P_tetraurelia:1[&&NHX:Ev=S:S=18:ND=18]):1[&&NHX:Ev=S:S=16:ND=16]):1[&&NHX:Ev=S:S=12:ND=12]):1[&&NHX:Ev=S:S=9:ND=9],P_tredecaurelia:1[&&NHX:Ev=S:S=10:ND=10]):1[&&NHX:Ev=S:S=6:ND=6]):1[&&NHX:Ev=S:S=4:ND=4]):1[&&NHX:Ev=S:S=1:ND=1],T_thermophila:1[&&NHX:Ev=S:S=2:ND=2])[&&NHX:Ev=S:S=0:ND=0];')

colorNodes(wgd1, 0)
ts = TreeStyle()
#ts.show_leaf_name = False
ts.show_scale = False
#wgd1.show(tree_style = ts)

st = basest.copy()
# read table with branch attributes
atrF = '/home/dsellis/data/IES/analysis/tables/ratesPerBranch.dat'
bratr = {}
f = open(atrF, 'r')
f.readline() # header
for line in f:
    line = line.rstrip()
    (fromP, toP, over, under) = line.split("\t")
    bratr[(fromP, toP)] = {'over': over, 'under': under}

colorNodes(st, 1)
# plot species tree
st.render('/home/dsellis/data/IES/analysis/figures/sptree.png', tree_style = ts)

glst = basest.copy()
colorNodes(glst, 0)
addattr(glst, bratr)
#glst.show(tree_style = ts)
glst.render('/home/dsellis/data/IES/analysis/figures/sptreeGL.png', tree_style = ts)
#quit()
agest = basest.copy()
for node in agest.traverse():
    figF = "/home/dsellis/data/IES/analysis/figures/length." + node.S + ".png"
    if os.path.isfile(figF):
        node.add_face(ImgFace(figF), column = 0)

# species tree with length plots
#agest.show(tree_style = ts)
#quit()
# Gene family trees #
#####################

geneFamily = sys.argv[1]

#  load gene family tree
treeF = '/home/dsellis/data/IES/analysis/phyldogT1/results/' + geneFamily + '.ReconciledTree'
t = Tree(treeF)
if not os.path.isfile(treeF):
    quit("charMatsF is not a file")

# # set styles
t = colorNodes(t, 1)

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

# load mean ancestral states
# asrF = '/home/dsellis/data/IES/analysis/tables/avNodeProb.dat';
# f = open (asrF, 'r')
# header = f.readline()

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
                    iesmotif.append( [int(beginMSA), int(endMSA), "[]", None, 10, "black", "red", "arial|8|black|" + iesId] )
    if iesmotif:
        seqFace = SeqMotifFace(seq = None, motifs = iesmotif, gap_format = "line")
        l.add_face(seqFace, 0, "aligned")
#           print("  " + geneId)
# # load alignment add biopython
# #alnF = '/home/dsellis/data/analysis/msas/filtered/cluster.10000.nucl.fa'

t.show()
"""

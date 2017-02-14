#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree, SeqMotifFace, ImgFace, SeqGroup, BarChartFace, PieChartFace, CircleFace
from pyies.userOptions import basePath
import os.path
import sys, getopt
import string
from collections import defaultdict
from pyies.functions import *

# plot gene family phylogenies with following styles:
#   gene tree with PHYLDOG species numbering and duplication nodes annotated
#   gene tree with IES ancestral state probabilities and observed presence/absence
#   gene tree with nucleotide alignment and schematic of IESs
# parameters

includedGeneFamilies = []
includeGF = ""
charMatFile = "analysis/iesdb/charMats.tab"
outputPath = ""
plotStyle = '1'
analysis = '2'
phyldogResultsPath = "analysis/phyldogT" + analysis + "/results"
ancNodeProbFile = "analysis/tables/avNodeProb" + analysis + ".dat"
nodeDictionaryFile = "analysis/tables/nodeDictionary" + analysis + ".dat"
homiesLinkFile = "analysis/tables/homIES" + analysis + ".columns.link"

def drawTree(outputFile):
    """ Print or save output and end script."""
    if outputFile:
        t.render(outputFile)
    else:
        t.show()

usage = """
usage:

./plotTree.py [OPTIONS]

    where OPTIONS can be any of the following:
    -g <geneFamily>
    -a [1|2|3] Choice of species tree analysis to use (Default: 2)
    -o <outputPath> Path to output the file named geneFamily.a.s.png
    -s [1|2|3] plot style
    -h this help screen
    -i file with gene families to include (ignoring option -g)

The coding of plot style used is the following:
1.  gene tree with PHYLDOG species numbering and duplication nodes annotated
2.  gene tree with IES ancestral state probabilities and observed presence/absence
3.  gene tree with nucleotide alignment and schematic of IESs

The coding for species tree is:
1. PHYLDOG species tree
2. concatenate with simple evolution model
3. concatenate with best codon models
"""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hg:a:o:i:s:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == "-g":
        includedGeneFamilies = [arg]
    elif opt == "-o":
        outputPath = arg
    elif opt == "-i":
        includeGF = arg
    elif opt == "-a":
        if arg in ['1','2','3']:
            analysis = arg
        else:
            print("Unknown value for option analysis (-a)")
            print(usage)
            sys.exit(1)
    elif opt == "-s":
        if arg in ['1','2','3']:
            plotStyle = arg
        else:
            print("Unknown value for option plotStyle (-s)")
            print(usage)
            sys.exit(1)

phyldogResultsPath = "analysis/phyldogT" + analysis + "/results"
ancNodeProbFile = "analysis/tables/avNodeProb" + analysis + ".dat"
nodeDictionaryFile = "analysis/tables/nodeDictionary" + analysis + ".dat"
homiesLinkFile = "analysis/tables/homIES" + analysis + ".columns.link"

# if parameter defined load a list of gene families to include from the analysis
if includeGF:
    with open(includeGF, 'r') as f:
        for line in f:
            includedGeneFamilies = [line.rstrip() for line in f]

# create dictionaries phyldog for translation:
rb2ph = {} # revBayes to phyldog node notation
hiesL = {} # homologous IES id to column number

# load nodeDictionary
nd = open(os.path.join(basePath, nodeDictionaryFile), 'r')
header = nd.readline()
for line in nd:
    line = line.rstrip()
    (cluster, r, phyldog, rb) = line.split("\t")
    if cluster in includedGeneFamilies:
        rb2ph[(cluster, rb)] = phyldog

# load link file connecting iesColumnId with column numbering in revBayes
lf = open(os.path.join(basePath, homiesLinkFile), 'r')
header = lf.readline()
for line in lf:
    line = line.rstrip()
    (gf, homIES, column) = line.split("\t")
    if gf in includedGeneFamilies:
        hiesL[(gf, homIES)] = int(column)

# load character matrices
f = open (os.path.join(basePath, charMatFile), 'r')
header = f.readline()
gfGenes  = defaultdict(list) # dictionary key is gene family value list of genes
gfhomIES = defaultdict(set) # dictionary key is gene family value list of homologous IESs
charMat = {}
for line in f:
    line = line.rstrip()
    (cluster, column, geneId, begin, end, ies, iesId, beginMSA, endMSA) = line.split("\t")
    gfGenes[cluster].append(geneId)
    gfhomIES[cluster].add(column)
    charMat[(cluster, column, geneId)] = [begin, end, ies, iesId, beginMSA, endMSA]

# for each gene family
for geneFamily in includedGeneFamilies:
    print("Ploting gene family: " + geneFamily)
    #  load gene family tree
    outputFile = os.path.join(outputPath, geneFamily + '.' + analysis + '.' + plotStyle + '.png')
    geneFamFile = geneFamily + ".ReconciledTree"
    treeF = os.path.join(basePath, phyldogResultsPath, geneFamFile)

    t = Tree(treeF)

    if plotStyle == '1':
        t = colorNodes(t, 1)
        drawTree(outputFile)
    elif plotStyle == '2':
        # load mean ancestral states
        # avNodeProb has homIES not with ids but with increment numbering
        fa = open(os.path.join(basePath, ancNodeProbFile), 'r')
        header = fa.readline()
        anc = defaultdict(dict) # defaultdictionary with key revbayes node id and value list of prob. presence
        for line in fa:
            line = line.rstrip()
            (cluster, node, iesColumn, presence) = line.split("\t")
            if cluster == geneFamily:
                anc[rb2ph[(cluster, node)]][iesColumn] = float(presence)

        # annotate tree
        for node in t.traverse():
            for homIES, presence in anc[node.ND].items():
                p = 100*presence
                a = 100 - p
                pf = PieChartFace([p, a], 10, 10, ["black", "silver"])
                if node.is_leaf():
                    pass
                else:
                    column = int(homIES) + 1
                    node.add_face(pf, column, "float")

        for leaf in t:
            geneId = leaf.name
            for homIES in gfhomIES[geneFamily]:
                cf = CircleFace(10, "white")
                (begin, end, ies, iesId, beginMSA, endMSA) = charMat[(geneFamily, homIES, geneId)]
                if ies == '?':
                    cf = CircleFace(10, "silver", label = "?")
                elif ies == '1':
                    cf = CircleFace(10, "black")
                elif ies == '0':
                    cf = CircleFace(10, "LightGrey")
                else:
                    sys.exit(1)
                column = hiesL[(geneFamily, homIES)] + 1
                leaf.add_face(cf, column, "aligned")
        drawTree(outputFile)
    elif plotStyle == '3': # plot with MSA
        # load nucleotide sequences for all genes!
        nuclAlnFile = os.path.join(basePath, 'analysis', 'msas', 'filtered', 'cluster.' + geneFamily + '.nucl.fa')
        seqs = SeqGroup(sequences = nuclAlnFile, format = "fasta")

        for leaf in t:
            geneId = leaf.name
            seq = seqs.get_seq(geneId)
            seq = seq.translate(None, string.ascii_lowercase) # keep only CDS
            iesmotif = [[1, len(seq), "line", 2, 5, None, None, None]]
            for homIES in gfhomIES[geneFamily]:
                (begin, end, ies, iesId, beginMSA, endMSA) = charMat[(geneFamily, homIES, geneId)]
                if ies == '?':
                    if beginMSA == 'NA':
                        iesmotif.append([int(begin), int(end),"()", 10, 10, "red", "black", "arial|8|black|?"])
                    else:
                        iesmotif.append([int(begin), int(end),"()", 10, 10, "red", "black", "arial|8|black|?"])
                elif ies == '1':
                    iesmotif.append([int(beginMSA), int(endMSA),"[]", 10, 10, "black", "red", "arial|8|black|" + iesId])
                elif ies == '0':
                    iesmotif.append([int(begin), int(end), "[]", 10, 10, "silver", "silver", None])
                else:
                    quit(1)
            seqFace = SeqMotifFace(seq = seq, motifs = iesmotif, gap_format = "blank", seq_format = "line")
            leaf.add_face(seqFace, 0, "aligned")
        drawTree(outputFile)


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

t.show()"""

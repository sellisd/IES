#!/usr/bin/python
from __future__ import print_function
import re
from ete3 import Tree, NodeStyle, TreeStyle, TextFace
from pyies.functions import readPalette9
import sys, getopt


inputfile = ''
outputfile = ''
usage = "./plotIESPhylogeny.py -i <inputfile> -o <outputfile> -t title"
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:t:",["ifile=","ofile=","cutoff="])
except getopt.GetoptError:
    print(usage)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt in ("-i", "--ifile"):
        inputfile = arg
    elif opt in ("-o", "--ofile"):
        outputfile = arg
    elif opt in ("-t", "--title"):
        title = arg


t = Tree(inputfile)
#'/home/dsellis/data/IES/analysis/mies/tempfiles/cl35886.aln.contree')
cp = readPalette9()

#cp  = ('1','2','3','4','5','6','7','8','9','10')
for node in t:
    m = re.search('(ppr|pbi|pte|ppe|pse|poc|ptr|pso|pca|tth)', node.name)
    if(m):
        nstyle = NodeStyle()
        nstyle["size"] = 10
        abr = m.group(0)
        if abr == 'ppr':
            nstyle['bgcolor'] = cp[0]
        elif abr == 'pbi':
            nstyle['bgcolor'] = cp[1]
        elif abr == 'pte':
            nstyle['bgcolor'] = cp[2]
        elif abr == 'ppe':
            nstyle['bgcolor'] = cp[3]
        elif abr == 'pse':
            nstyle['bgcolor'] = cp[4]
        elif abr == 'poc':
            nstyle['bgcolor'] = cp[5]
        elif abr == 'ptr':
            nstyle['bgcolor'] = cp[6]
        elif abr == 'pso':
            nstyle['bgcolor'] = cp[7]
        elif abr == 'pca':
            nstyle['bgcolor'] = cp[8]
        elif abr == 'tth':
            nstyle['bgcolor'] = cp[9]
        else:
            nstyle['bgcolor'] = "#000000"
        node.set_style(nstyle)


ts = TreeStyle()
#ts.show_leaf_name = False
ts.mode = 'c'
ts.title.add_face(TextFace(title, fsize=20), column=0)
#t.show(tree_style = ts)
t.render(outputfile, tree_style = ts)

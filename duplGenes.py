#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree, NodeStyle, SeqMotifFace, ImgFace, TreeStyle
from pyies.functions import placeDupl, isDescendant, choseAnc, up2sp, choseAncList
from os import listdir
import os.path
from collections import Counter

#Count the per branch number of duplication events.

spt = Tree('((P_caudatum:1[&&NHX:Ev=S:S=3:ND=3],((P_sexaurelia:1[&&NHX:Ev=S:S=7:ND=7],P_sonneborni:1[&&NHX:Ev=S:S=8:ND=8]):1[&&NHX:Ev=S:S=5:ND=5],(((P_pentaurelia:1[&&NHX:Ev=S:S=13:ND=13],P_primaurelia:1[&&NHX:Ev=S:S=14:ND=14]):1[&&NHX:Ev=S:S=11:ND=11],(P_biaurelia:1[&&NHX:Ev=S:S=15:ND=15],(P_octaurelia:1[&&NHX:Ev=S:S=17:ND=17],P_tetraurelia:1[&&NHX:Ev=S:S=18:ND=18]):1[&&NHX:Ev=S:S=16:ND=16]):1[&&NHX:Ev=S:S=12:ND=12]):1[&&NHX:Ev=S:S=9:ND=9],P_tredecaurelia:1[&&NHX:Ev=S:S=10:ND=10]):1[&&NHX:Ev=S:S=6:ND=6]):1[&&NHX:Ev=S:S=4:ND=4]):1[&&NHX:Ev=S:S=1:ND=1],T_thermophila:1[&&NHX:Ev=S:S=2:ND=2])[&&NHX:Ev=S:S=0:ND=0];')

baseP = '/home/dsellis/data/IES/analysis/phyldog/results/'
treeFs = [i for i in listdir(baseP) if i[-15:] == '.ReconciledTree']
print("from\tto\tduplications")
dup = Counter()
for i in treeFs[1:100]:
    geneFamily = i[0:-15]
    treeF = os.path.join(baseP, i)
    t = Tree(treeF)
    if not os.path.isfile(treeF):
        quit("treeF is not a file")
    dup = dup + placeDupl(t, spt)
for i in dup:
    print(str(i[0]) + "\t" + str(i[1]) + "\t" + str(dup[i]))

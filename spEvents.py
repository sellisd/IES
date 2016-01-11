#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree

"""
find which phyldog nodes correspond to which speciation event
"""

phyldogPath = '/home/dsellis/data/IES_data/msas/phyldog/results/'
inputF = open('/home/dsellis/data/IES_data/msas/asr/geneFamilies.dat', 'r')
clusters = inputF.readlines()
clusters = [i.rstrip() for i in clusters]
print('cluster\tnodeP\tspEvent')

for cluster in clusters:
    fileNameString = phyldogPath+str(cluster)+'.ReconciledTree'
    #    print(fileNameString)
    t = Tree(fileNameString)
    for node in t:
        if(node.Ev == 'S'): # if speciation event
            print(cluster + "\t" + node.ND + "\t" + node.S)


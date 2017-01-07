#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree
import sys, getopt
import os.path

# Find which phyldog nodes correspond to which speciation event.

phyldogPath = ''
geneFamilies = ''
usage = "./spEvents.py -p <phyldogResultsPath> -g <geneFamilies>"
try:
    opts, args = getopt.getopt(sys.argv[1:],"hp:g:")
except getopt.GetoptError:
    print(usage)
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit(1)
    elif opt == '-p':
        phyldogPath = arg
    elif opt == '-g':
        geneFamilies = arg

inputF = open(geneFamilies, 'r')
clusters = inputF.readlines()
clusters = [i.rstrip() for i in clusters]
print('cluster\tnodeP\tspEvent')

for cluster in clusters:
    fileNameString = os.path.join(phyldogPath, str(cluster)+'.ReconciledTree')
    t = Tree(fileNameString)
    for node in t.traverse():
        if(node.Ev == 'S'): # if speciation event
            print(cluster + "\t" + node.ND + "\t" + node.S)

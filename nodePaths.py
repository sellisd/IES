#!/usr/bin/python
from __future__ import print_function
from ete3 import Tree

"""
read a gene tree (nhx format) and for each pair of speciation nodes print a path of ancestor to descendant.
"""

# Functions
def path2anc(nodeO):
    """Find path to root
    
    PARAMETERS: reference to a node object
    RETURN    : a list of node objects on the path from the root to the input node
    """
    path = []
    while nodeO.up:
        path.append(nodeO)
        nodeO = nodeO.up
    path.append(nodeO) # add root
    return path

def subpaths(pairS, L):
    """
    Find subpaths in node list L
    Given pairs of speciation events and a list of nodes returns subpaths starting and ending on given pairs
    """
    subpaths = []
    for p in pairS:
        subpath = []
        inSubpath = False
        for n in L:
            if inSubpath:
                subpath.append(n)
                if n.Ev == 'S' and str(p[1]) == n.S:
                    inSubpath = False
                    break
            else:
                if n.Ev == 'S' and str(p[0]) == n.S:
                    subpath.append(n)
                    inSubpath = True
        else:
            subpath = [] # discard if not completed
        if subpath:
            subpaths.append(subpath)
    return subpaths

#pairs of speciation nodes

spNodePairsF = '/home/dsellis/data/IES/analysis/tables/spNodePairs.dat'
pairS = [i.rstrip().split(" ") for i in open(spNodePairsF, 'r')]
#pairS = ((0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(2,3),(2,4),(2,5),(2,6),(4,5),(4,6))
phyldogPath = '/home/dsellis/data/IES/analysis/phyldog/results/'

inputF = open('/home/dsellis/data/IES/analysis/asr/geneFamilies.dat', 'r')
clusters = inputF.readlines()
clusters = [i.rstrip() for i in clusters]
print('cluster\tfrom\tto\tpath')

for cluster in clusters:
    fileNameString = phyldogPath+str(cluster)+'.ReconciledTree'
    # print(fileNameString)
    t = Tree(fileNameString)
    # print(t.get_ascii(attributes=['ND','name']))
    for leaf in t:
        L = path2anc(leaf)
        L.reverse();
        # print(list([n.ND for n in L]))
        # print(list([n.Ev for n in L]))
        # print(list([n.S for n in L]))
        Sbs = subpaths(pairS,L)
        for i in Sbs:
            toP = list([j.ND for j in i])
            fromS = i[0].S
            toS = i[-1].S
            print('\t'.join([str(cluster), fromS,toS, ','.join(toP)]))




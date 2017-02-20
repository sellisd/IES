#!/usr/bin/python
from __future__ import print_function
from pyies.functions import makeNodeTable
from pyies.userOptions import basePath
import os.path

# read list of gene families
rbNodeIndexes = os.path.join(basePath, 'analysis', 'asr1', 'rbNodeIndexes')
phyldogResults = os.path.join(basePath, 'analysis', 'phyldogT1', 'results')

with open(os.path.join(basePath, 'analysis', 'iesdb', 'geneFamilydb.dat'), 'r') as f:
#f = open('workfile', 'r')
    #with open('dd', 'r') as f:
    for line in f:
        line = line.rstrip()
        (gfId, seqNo) = line.split("\t")[0:2]
        if seqNo == 'NA':
            continue
        rbTree = os.path.join(rbNodeIndexes, 'nodeIndex.' + str(gfId) + '.tre')
        phTree = os.path.join(phyldogResults, str(gfId) + '.ReconciledTree')
        if(os.path.isfile(rbTree) & os.path.isfile(phTree)):
             L = makeNodeTable(phyldogTreeF = phTree, revBayesTreeF = rbTree)
             for (rb, ph) in L:
                 print("\t".join([str(gfId), rb, ph]))

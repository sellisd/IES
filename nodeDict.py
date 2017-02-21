#!/usr/bin/python
from __future__ import print_function
from pyies.functions import makeNodeTable
from pyies.userOptions import basePath
import os.path

"""Create dictionary of nodeIds."""

print("Reading gene families.")
geneFamilies = []
with open(os.path.join(basePath, 'analysis', 'iesdb', 'geneFamilydb.dat'), 'r') as f:
    for line in f:
        line = line.rstrip()
        (gfId, seqNo) = line.split("\t")[0:2]
        if seqNo == 'NA':
            continue
        else:
            geneFamilies.append(gfId)

for asrRun in (1,2,3):
    print("Processing run " + str(asrRun))
    rbNodeIndexes = os.path.join(basePath, 'analysis', 'asr' + str(asrRun), 'rbNodeIndexes')
    phyldogResults = os.path.join(basePath, 'analysis', 'phyldogT' + str(asrRun), 'results')
    out = os.path.join(basePath, 'analysis', 'tables', 'nodeDictionary' + str(asrRun) + '.dat')
    with open(out, 'w') as fout:
        fout.write("\t".join(['cluster', 'r', 'phyldog', 'rb\n']))
        for gfId in geneFamilies:
            rbTree = os.path.join(rbNodeIndexes, 'nodeIndex.' + str(gfId) + '.tre')
            phTree = os.path.join(phyldogResults, str(gfId) + '.ReconciledTree')
            if(os.path.isfile(rbTree) & os.path.isfile(phTree)):
                L = makeNodeTable(phyldogTreeF = phTree, revBayesTreeF = rbTree)
                for (rb, ph) in L:
                    fout.write("\t".join([str(gfId), 'NA', ph, rb + '\n']))

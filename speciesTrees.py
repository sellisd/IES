#!/usr/bin/env python
from __future__ import print_function
import os.path
from ete3 import TreeStyle
from pyies.userOptions import basePath
from pyies.functions import phyldogSpeciesTree

#generate species trees
outgroupName = "Tetrahymena_thermophila"

brlenT1 = os.path.join(basePath, 'analysis', 'sgf', 'topoConstrSimple.treefile')
brlenT2 = os.path.join(basePath, 'analysis', 'sgf', 'concatSimple.nexus.treefile')
brlenT3 = os.path.join(basePath, 'analysis', 'sgf', 'concat.nexus.treefile')
brlenT3b = os.path.join(basePath, 'analysis', 'sgf', 'concatBrlen.treefile')
phT1 = os.path.join(basePath, 'analysis', 'phyldogT1', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree')
phT2 = os.path.join(basePath, 'analysis', 'phyldogT2', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree')
phT3 = os.path.join(basePath, 'analysis', 'phyldogT3', 'results', 'OutputSpeciesTree_ConsensusNumbered.tree')

t1 = phyldogSpeciesTree(phT1, brlenT1, outgroupName)
t2 = phyldogSpeciesTree(phT2, brlenT2, outgroupName)
t3 = phyldogSpeciesTree(phT3, brlenT3, outgroupName, codon = True)
t3b = phyldogSpeciesTree(phT3, brlenT3b, outgroupName)

t1.write(outfile = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree1.nhx'), features=['ND'])
t2.write(outfile = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree2.nhx'), features=['ND'])
t3.write(outfile = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree3.nhx'), features=['ND'])
t3b.write(outfile = os.path.join(basePath, 'analysis', 'iesdb', 'speciesTree3b.nhx'), features=['ND'])

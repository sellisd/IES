#!/usr/bin/env python
from __future__ import print_function
import os.path
from ete3 import TreeStyle
from pyies.userOptions import basePath
from pyies.functions import phyldogSpeciesTree

#generate species trees
outpgroupName = "Tetrahymena_thermophila"

brlenT1 = os.path.join(basePath, 'analysis', 'sgf', '')
brlenT2 = os.path.join(basePath, 'analysis', 'sgf', '')
brlenT3 = os.path.join(basePath, 'analysis', 'sgf', '')

phT1 = os.path.join(basePath, 'analysis', 'phyldogT1', 'results')
phT2 = os.path.join(basePath, 'analysis', 'phyldogT2', 'results')
phT3 = os.path.join(basePath, 'analysis', 'phyldogT3', 'results')

t1 = phyldogSpeciesTree(phT1, brlenT1, outgroupName)
t2 = phyldogSpeciesTree(phT2, brlenT2, outgroupName)
t3 = phyldogSpeciesTree(phT3, brlenT3, outgroupName)

t1.write(outfile="speciesTree1.nhx")
t2.write(outfile="speciesTree2.nhx")
t3.write(outfile="speciesTree3.nhx")

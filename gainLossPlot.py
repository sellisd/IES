#!/usr/bin/evn python
from __future__ import print_function
from __future__ import division
from pyies.userOptions import basePath
from ete3 import TextFace, TreeStyle, NodeStyle
import pandas as pd
outputF = ""
inputGLF = os.path.join(basePath, 'analysis', 'tables', 'gainLossSum1.dat')
inputTF = "/Volumes/WDC/data/IES/analysis/sgf/topoConstrSimple.treefile"
gl = pd.read_csv(inputGLF)
t = Tree(inputTF)

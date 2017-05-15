#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import pandas as pd
from pyies.userOptions import basePath
import os.path
import glob
from ete3 import Tree
from pyies.functions import NDS

for tr in [1,2,3]:
  fp = os.path.join(basePath, "analysis", "phyldogT" + str(tr), "results", "*.ReconciledTree")
  for f in glob.glob(fp):
    print(NDS(Tree(f)))
    quit()

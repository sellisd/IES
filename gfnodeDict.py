#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import pandas as pd
from pyies.userOptions import basePath
import os.path
import glob
import re
from ete3 import Tree
from pyies.functions import NDS

for tr in [1,2,3]:
  with open(os.path.join(basePath, "analysis", "tables", "gfnodeDict" + str(tr) + ".tsv"), 'w') as fout:
      fp = os.path.join(basePath, "analysis", "phyldogT" + str(tr), "results", "*.ReconciledTree")
      fout.write("\t".join(['geneFamily', 'ND', 'S\n']))
      for f in glob.glob(fp):
        match = re.search(r"^.*?(\d+)\.ReconciledTree$", f)
        geneFamily = match.group(1)
        print(str(tr) + "/" + geneFamily + "\r")
        for (k,v) in NDS(Tree(f)).items():
            fout.write("\t".join([geneFamily, k, v + "\n"]))

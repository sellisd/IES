#!/usr/bin/python
from __future__ import print_function
from pyies import functions
import sys, os, fnmatch, re

# parse multiple PHYLDOG output and compare gene tree likelihoods
runPaths = sys.argv[1:]
loglk = {}
runs = []
for path in runPaths:
    files = []
    m = re.search('/phyldogT(\d+)/', path)
    if m:
        run = 'run' + m.group(1)
    else:
        quit("pathname error:" + path)
    for f in os.listdir(path):
        if fnmatch.fnmatch(f, "Client_[12].out"):
            files.append(open(os.path.join(path,f), 'r'))
    loglk = functions.parseClientOut(files, loglk, run)
    runs.append(run)


print("geneFamily", "\t".join(runs), sep = "\t", end = "\n")
geneFamilies = set([g for (g, r) in loglk.keys()])

for gf in geneFamilies:
    print(gf, end = "\t")
    for r in runs:
        print(loglk[(gf, r)], end = "\t")
    print()



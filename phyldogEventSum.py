#!/usr/bin/python
from __future__ import print_function
from collections import Counter
import sys
import re

# Calculate number of events per species tree branch.

dupl = Counter() #./phyldogEventSum.py pathToPhyldog/results/*_Events.txt
loss = Counter()
nodes = set()
for fn in sys.argv[1:]:
    f = open(fn, 'r')
    for line in f:
        line = line.rstrip()
        m = re.search('event\((\d+).*/(\d+)\.opt",(duplication|loss)\)', line)
        if(m):
            (node, geneFamily, event) = m.group(1,2,3)
            node = int(node)
            if(event == 'duplication'):
                dupl[node] += 1
                nodes.add(node)
            else:
                loss[node] += 1
                nodes.add(node)
        else:
            print(line)

print("speciesNode\tduplications\tlosses")
for d in sorted(nodes):
    print(d, dupl[d], loss[d], sep = "\t")

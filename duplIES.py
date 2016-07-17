#!/usr/bin/python
from __future__ import print_function
import sys
from collections import Counter
from string import maketrans

"""Find identical IES within one species."""

def revcomp(a):
	r = a[::-1] # reverse
	translation = maketrans(("ATCGatcg"),("TAGCtagc"))
	rc = r.translate(translation)
	return(rc)

file = open(sys.argv[1], "r")

file.readline() #header
dupl = Counter()
for line in file:
	line.rstrip()
	(IESid, scaffold, altSeqNo, startLocs, endLocs, upstreamFlank, downstreamFlank, length, isFloating, merged, inCDS, inInter, inIntron, front, back, seq) = line.split()
	if isFloating == "0" and merged == "0":
		if seq in dupl or revcomp(seq) in dupl:
			dupl[seq] += 1
		else:
			dupl[seq] = 0

#loop through dictionary and find how many with more than 2 identical copies
for i in dupl.keys():
	if dupl[i] > 1:
		print(str(dupl[i]) + ' ' + i)
#!/usr/bin/python
from __future__ import print_function
import sys

# Modify option files for PHYLDOG.

for fn in sys.argv[1:]: # ./addPhylOpt.py PATH2OPT/*.opt
    f = open(fn,'a')
    f.write("\noutput.events.file=$(RESULT)$(DATA)_Events.txt\n")

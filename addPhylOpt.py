#!/usr/bin/python
from __future__ import print_function
import sys
""" Modify option files for PHYLDOG."""
for fn in sys.argv[1:]:
    f = open(fn,'a')
    f.write("\noutput.events.file=$(PATH)$(DATA)_Events.txt\n")
    

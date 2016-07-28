# sim_max.py
#
# Find the greatest simulated capacities.
#
# This code is experimental, and error-handling is primitive.
#
# Copyright 2013, NICTA.  See COPYRIGHT for license details.

#!/usr/bin/env python

import sys

ymax= -1

for l in sys.stdin:
    bits= l.strip().split(' ')
    if len(bits) == 0: continue
    if bits[0] == '#': continue

    x, y, ye= map(float, bits)

    if y+ye > ymax: ymax= y+ye

print "%.2e %.12e" % (ymax, ymax)

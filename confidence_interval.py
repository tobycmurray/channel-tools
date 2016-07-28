# confidence_interval.py
#
# Find the estimated confidence interval for a given capacity, using the
# output of sample_error.
#
# This code is experimental, and error-handling is primitive.
#
# Copyright 2013, NICTA.  See COPYRIGHT for license details.

#!/usr/bin/env python

import sys

meas_low=  float(sys.argv[1])
meas_high= float(sys.argv[2])

real_caps= []
bounds= dict()

# Find the min and max observed capacities for each real capacity.
for l in sys.stdin:
    bits= l.strip().split(' ')
    if len(bits) == 0: continue
    if bits[0] == '#': continue

    creal= float(bits[0])
    cmeas= float(bits[1])

    if not creal in bounds:
        real_caps.append(creal)
        bounds[creal]= [cmeas, cmeas]

    if cmeas < bounds[creal][0]:
        bounds[creal][0]= cmeas
    if cmeas > bounds[creal][1]:
        bounds[creal][1]= cmeas

# Find the first interval that contains the lowest measured capacity
ilow= None
for i in xrange(len(real_caps)):
    c= real_caps[i]
    if meas_low <= bounds[c][1]:
        ilow= i
        break

if ilow is None:
    print "Out of range (high)"
    sys.exit(1)
else:
    if ilow > 0:
        ca= real_caps[ilow]
        cb= real_caps[ilow-1]
        x= (bounds[ca][1] - meas_low) / \
           (bounds[ca][1] - bounds[cb][1])
        clow= x * cb + (1-x) * ca
    else:
        c= real_caps[ilow]
        
# Find the last interval that contains the highest measured capacity
ihigh= None
for i in reversed(xrange(len(real_caps))):
    c= real_caps[i]
    if bounds[c][0] <= meas_high:
        ihigh= i
        break

if ihigh is None:
    print "Out of range (low)"
    sys.exit(1)
else:
    if ilow < len(real_caps) - 1:
        ca= real_caps[ilow]
        cb= real_caps[ilow+1]
        x= (meas_high - bounds[ca][0]) / \
           (bounds[cb][0] - bounds[ca][0])
        chigh= x * cb + (1-x) * ca
    else:
        chigh= real_caps[ihigh]

print (clow+chigh)/2, clow, chigh

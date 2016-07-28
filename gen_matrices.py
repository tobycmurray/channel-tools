# gen_matrices.py
#
# Generate all channel matrices for the given results hierarchy.
#
# This code is experimental, and error-handling is primitive.
#
# Copyright 2013, NICTA.  See COPYRIGHT for license details.

#!/usr/bin/env python

import os.path
from subprocess import Popen, PIPE
import sys

dataroot= sys.argv[1]
scriptroot= os.path.dirname(sys.argv[0])

filter_samples= os.path.join(scriptroot, "filter_samples")
channel_matrix= os.path.join(scriptroot, "channel_matrix")

pipes= []

stats= file(os.path.join(dataroot, "analysis", "stats"))
for l in stats:
    bits= l.strip().split(' ')
    chip, chan, cm, ts, cmin, cmax, rmin, rmax= bits[0:8]
    min_count= bits[-1]

    output= os.path.join(dataroot, "matrices", \
        "%s.%s.%s.%s.cm" % (chip, chan, cm, ts))

    print "Generating %s\n" % output

    cmdline= 'find %s/%s/%s/%s/TS_%s -name "*.xz" | ' % \
               (dataroot, chip, chan, cm, ts ) + \
             'xargs xzcat | ' + \
             '%s %s %s %s %s | ' % \
               (filter_samples, cmin, cmax, rmin, rmax) + \
             '%s %s %s %s %s > /dev/null 2>&1' % \
               (channel_matrix, output, cmin, cmax, min_count)

    p= Popen(cmdline, shell=True, stdout=PIPE)
    pipes.append(p)

for p in pipes:
    p.communicate()

stats.close()

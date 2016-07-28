# deficit.py
#
# Output the number of samples required for each modulation to hit the target.
#
# This code is experimental, and error-handling is primitive.
#
# Copyright 2013, NICTA.  See COPYRIGHT for license details.

#!/usr/bin/env python

import sys

target= int(sys.argv[1])

for l in sys.stdin:
    bits= l.strip().split(' ')
    if len(bits) == 2:
        c, count= map(int, bits)
        if count < target:
            print c, target - count

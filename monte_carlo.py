# monte_carlo.py
#
# Run error simulations for the given channel.
#
# This code is experimental, and error-handling is primitive.
#
# Copyright 2013, NICTA.  See COPYRIGHT for license details.

#!/usr/bin/env python

import os.path
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 10:
    print >>sys.stderr, "Usage: %s <data root> <chip> <channel> " \
        "<countermeasure> <timeslice> <divide limit> <size limit> " \
        "<precision> <runs>" % \
        sys.argv[0]
    sys.exit(1)

jobs= 0
cpuinfo= file("/proc/cpuinfo", "r")
for l in cpuinfo:
    if re.match("^processor", l):
        jobs+= 1
if jobs < 1:
    jobs= 1

rootdir= os.path.dirname(sys.argv[0])

# Paths to utilities
channel_hist= os.path.join(rootdir, "channel_hist")
channel_matrix= os.path.join(rootdir, "channel_matrix")
capacity= os.path.join(rootdir, "capacity")
sample_error= os.path.join(rootdir, "sample_error")
filter_samples= os.path.join(rootdir, "filter_samples")

dataroot, chip, channel, countermeasure= sys.argv[1:5]
timeslice, divlimit, limit= map(int, sys.argv[5:8]) 
epsilon= float(sys.argv[8])

runs= int(sys.argv[9])

print "Monte-carlo simulation of hypothetical empty-channel bandwidth."
print "Chip: %s, Channel: %s, Countermeasure: %s, Timeslice: %d" % \
    (chip, channel, countermeasure, timeslice)
print "Error limit: %.2e, Simulating %d matrices per size" % \
    (epsilon, runs)
print "Experimental data in %s" % os.path.abspath(dataroot)

run_name= "%s.%s.%s.%d" % (chip, channel, countermeasure, timeslice)

divisors= []
d= 1
for i in xrange(divlimit):
    divisors.append(d)
    d*= 2
sizes= [limit / d for d in divisors]

print "Sampling at sizes [" + ", ".join(map(str,sizes)) + "]"

channel_dir= os.path.join(dataroot, chip, channel)
cm_dir= os.path.join(channel_dir, countermeasure)
runs_dir= os.path.join(cm_dir, "TS_%d" % timeslice)

range_str= None
info= file(os.path.join(channel_dir, "info"), "r")
for l in info:
    m= re.match("^modulation range:\s*(\d+)\s*-\s*(\d+)", l)
    if m:
        range_str= m.group(1, 2)
del info
mod_range= (int(range_str[0]), int(range_str[1]))

limits= file(os.path.join(cm_dir, "limits"), "r")
for l in limits:
    m= re.match("^result range:\s*(\d+)\s*-\s*(\d+)", l)
    if m:
        range_str= m.group(1, 2)
del limits
res_range= (int(range_str[0]), int(range_str[1]))

print "Decompressing and collating samples...",
sys.stdout.flush()
samples_file= tempfile.mkstemp()
subprocess.call("find %s -name \"*.xz\" | xargs xzcat | " \
                "%s %s %s %s %s > %s" % \
    (runs_dir, filter_samples, mod_range[0], mod_range[1], \
     res_range[0], res_range[1], samples_file[1]), shell=True)
print "done."

matrix= tempfile.mkstemp()
subprocess.call("%s %s %d %d < %s" \
    % (channel_matrix, matrix[1], mod_range[0], mod_range[1], \
       samples_file[1]), shell=True)

print "Simulating noisy channel matrices"

sim_output= file(run_name + ".sim", "w")

for s in sizes:
    print "  Using %d samples" % s

    pipe= subprocess.Popen("%s %s %d %d %e -q" % (sample_error, \
        matrix[1], s, runs, epsilon), stdout=subprocess.PIPE,
        shell=True)

    sodata, sedata= pipe.communicate()
    for datum in sodata.strip().split('\n'):
        sim_output.write("%d %s\n" % (s, datum.strip()))

del(sim_output)
    
os.unlink(matrix[1])
os.unlink(samples_file[1])

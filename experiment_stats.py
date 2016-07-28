# experiment_stats.py
#
# Calculate statistics for all experimental data found in the given hierarchy.
#
# This code is experimental, and error-handling is primitive.
#
# Copyright 2013, NICTA.  See COPYRIGHT for license details.

#!/usr/bin/env python

import gzip
import os.path
from numpy import zeros, array
import re
from subprocess import Popen, PIPE
import sys

rootdir= os.path.dirname(sys.argv[0])

machine= False

if len(sys.argv) < 2:
    print >>sys.stderr, "Usage: %s <data root> [-m]\n" % sys.argv[0]
    sys.exit()

if len(sys.argv) > 2 and sys.argv[2] == '-m':
    machine= True

jobs= 0
cpuinfo= file("/proc/cpuinfo", "r")
for l in cpuinfo:
    if re.match("^processor", l):
        jobs+= 1
if jobs < 1:
    jobs= 1

root= os.path.abspath(sys.argv[1])
for chip in os.listdir(root):
    if not os.path.isdir(os.path.join(root, chip)):
        continue
    if chip == "simulation" or chip == "analysis" or \
       chip == "quarantine" or chip == "alternate":
        continue

    if not machine:
        print "chip %s" % chip
    chip_root= os.path.join(root, chip)

    channels= os.listdir(chip_root)
    for c in channels:
        cd= os.path.join(chip_root, c)
        if not os.path.isdir(cd):
            continue
        if c == "manual":
            continue

        rge_str= None
        info= file(os.path.join(cd, "info"), "r")
        for l in info:
            m= re.match("^modulation range:\s*(\d+)\s*-\s*(\d+)", l)
            if m:
                rge_str= m.group(1, 2)
        del info
        rge= (int(rge_str[0]), int(rge_str[1]))

        if not machine:
            print "  channel %s" % c
            print "  range %d - %d" % rge

        countermeasures= os.listdir(cd)
        for cm in countermeasures:
            cmd= os.path.join(cd, cm)
            if not os.path.isdir(cmd):
                continue
            if not machine:
                print "    countermeasure %s" % cm

            rge_str= None
            limits= file(os.path.join(cmd, "limits"), "r")
            for l in limits:
                m= re.match("^result range:\s*(\d+)\s*-\s*(\d+)", l)
                if m:
                    rge_str= m.group(1, 2)
            del limits
            res_rge= (int(rge_str[0]), int(rge_str[1]))

            timeslices= os.listdir(cmd)

            for ts_string in timeslices:
                if ts_string == "limits": continue

                tsd= os.path.join(cmd, ts_string)
                m= re.match("TS_(.*)", ts_string)
                ts= int(m.group(1))
                if not machine:
                    print "      timeslice %d" % ts

                runs= 0
                out_of_range= 0
                malformed= 0
                counts= zeros(rge[1] - rge[0] + 1)
                total_count= 0
                rmin= None
                rmax= None
                cmin= None
                cmax= None

                run_paths= []
                for build in os.listdir(tsd):
                    build_dir= os.path.join(tsd, build)
                    for run in os.listdir(build_dir):
                        run_path= os.path.join(build_dir, run)
                        if not os.path.isfile(run_path):
                            continue
                        run_paths.append(run_path)
                        runs+= 1

                job_paths= [[] for j in xrange(jobs)]
                j= 0
                for rp in run_paths:
                    job_paths[j].append(rp)
                    j= (j + 1) % jobs

                name= "%s.%s.%s.%d" % (chip, c, cm, ts)
                pipes= [Popen("xzcat %s | %s/summarise %d %d %d %d %s %s" % \
                            (" ".join(jp), rootdir, rge[0], rge[1], \
                             res_rge[0], res_rge[1], \
                             "%s.ool" % name, "%s.mal" % name), \
                            shell=True, stdout=PIPE) for jp in job_paths
                            if len(jp) > 0]

                for p in pipes:
                    stdoutdata, stderrdata= p.communicate()
                    del(p)

                    output= map(int, stdoutdata.strip().split(" "))

                    total_count+= output[0]
                    out_of_range+= output[1]
                    malformed+= output[2]
                    if rmin is None or output[3] < rmin:
                        rmin= output[3]
                    if rmax is None or rmax < output[4]:
                        rmax= output[4]
                    if cmin is None or output[5] < cmin:
                        cmin= output[5]
                    if cmax is None or cmax < output[6]:
                        cmax= output[6]
                    counts+= array(output[7:])

                if machine:
                    print "%s %s %s %d %d %d %d %d %d %d" \
                          " %d %d %d %d %d %d %d" % \
                        (chip, c, cm, ts, rge[0], rge[1], \
                         res_rge[0], res_rge[1], runs, total_count, \
                         out_of_range, malformed, rmin, rmax, \
                         cmin, cmax, min(counts))
                else:
                    print "        %d runs" % runs
                    print "        %d counts" % total_count
                    print "        %d out of range" % out_of_range
                    print "        %d malformed" % malformed
                    print "        %d min row" % rmin
                    print "        %d max row" % rmax
                    print "        %d min col" % cmin
                    print "        %d max col" % cmax
                    print "        %d minimum count" % min(counts)

                hist= file("%s.hist" % name, "w")
                for i in xrange(len(counts)):
                    print >>hist, "%d %d\n" % (rge[0] + i, counts[i])
                del(hist)

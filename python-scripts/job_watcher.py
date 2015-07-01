#!/usr/bin/env python

"""This Python script reports on running jobs."""

import sys
import os
import glob
import time

# list of directories to monitor
dlist = [ \
         '/ptmp/jhalvers/400_13_slow_insert_INT32_3', \
         '/ptmp/jhalvers/rings100_100_INT32_2', \
         '/ptmp/jhalvers/sllod_M200N800_rings_G0.0000001', \
         '/ptmp/jhalvers/sllod_M200N800_rings_G1e-6.5', \
         '/ptmp/jhalvers/sllod_M200N800_rings_G0.000001' \
        ]

# maximum number of time steps for LAMMPS
max_lammps = 2**31

def file_age(atime, now):
  mts = (now - atime) / 60.0
  hrs = mts / 60.0
  dys = mts / (60.0 * 24.0)
  if(mts < 60.0):
    s = '%.1f minutes' % mts
    return s
  elif(hrs < 24.0):
    s = '%.1f hrs' % hrs
    return s
  else:
    s = '%.1f days' % dys
    return s

s = ''
summary = ''
for d in dlist:

  # cd to directory d
  os.chdir(d)

  # grab all dump files report small and largest
  dumpfiles = glob.glob('rings*') + glob.glob('sllod*')
  dumpfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), dumpfiles)
  dumpfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  s = s + '\n\n' + d + '\n\n  Number of dump files: ' + str(len(dumpfiles)) + '\n'
  if(len(dumpfiles) > 1):
    for f in (dumpfiles[0], dumpfiles[-2], dumpfiles[-1]):
      sz = '%.1f MB' % (os.path.getsize(f) / 1e6)
      age = pypolymer.file_age(os.path.getctime(f), time.time())
      s = s + '    ' + f + '  ' + sz + '  ' + age + '\n'
    max_timestep = float(dumpfiles[-1][dumpfiles[-1].rindex('.') + 1:])
    percent_done = '%.1f' % (100 * max_timestep / max_lammps)
    summary = summary + percent_done + '% ' + d + '\n'

  logfiles = glob.glob('log.*')
  logfiles.sort(lambda u, v: cmp(int(u[u.index('.') + 1:]), int(v[v.index('.') + 1:])))
  if(len(logfiles) > 1):
    s = s + '\n  Number of log files: ' + str(len(logfiles)) + '\n'
    for f in logfiles:
      sz = '%.1f KB' % (os.path.getsize(f) / 1e3)
      age = pypolymer.file_age(os.path.getctime(f), time.time())
      s = s + '    ' + f + '  ' + sz + '  ' + age + '\n'

  stressfiles = glob.glob('stress.dat.*')
  stressfiles.sort(lambda u, v: cmp(int(u[u.index('.') + 5:]), int(v[v.index('.') + 5:])))
  if(len(stressfiles) > 1):
    s = s + '\n  Number of stress files: ' + str(len(stressfiles)) + '\n'
    for f in stressfiles:
      sz = '%.1f KB' % (os.path.getsize(f) / 1e3)
      age = pypolymer.file_age(os.path.getctime(f), time.time())
      s = s + '    ' + f + '  ' + sz + '  ' + age + '\n'

  restartfiles = glob.glob('restart.*')
  restartfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  if(len(restartfiles) > 1):
    s = s + '\n  Number of restart files: ' + str(len(restartfiles)) + '\n'
    for f in restartfiles:
      sz = '%.1f MB' % (os.path.getsize(f) / 1e6)
      age = pypolymer.file_age(os.path.getctime(f), time.time())
      s = s + '    ' + f + '  ' + sz + '  ' + age + '\n'

s = 'RZG jobs for jhalvers\n\n' + summary + s

# send results by email
cmd = 'echo "' + s + '" | mail -s "job_watcher results on RZG" jhalvers@rzg.mpg.de'
estat = os.system(cmd)

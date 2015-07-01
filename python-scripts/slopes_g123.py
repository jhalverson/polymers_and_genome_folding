#!/usr/bin/env python

"""This Python script computes the slope of g123 data as a function
   of time. It may be used to determine crossover times between the
   various regimes. Linear or logarithm scales may be considered."""

import sys
import os
import math
import numpy as np
import pypolymer

# g123 data file
base = '/people/thnfs/homes/halvers/'
fin = base + 'research/ring_polymers/g123_lammps_espresso_rings_combined/800/rings800_msd_origins_full.dat'
fin = base + 'research/ring_polymers/g123_lammps_espresso_rings_combined/400/rings400_msd_origins_full.dat'
fin = base + 'research/ring_polymers/g123_lammps_espresso_rings_combined/200/rings200_msd_origins_full.dat'
fin = base + 'research/ring_polymers/g123_lammps_espresso_rings_combined/100/rings100_msd_origins_full.dat'
fin = base + 'research/ring_polymers/g123_lammps_espresso_rings_combined/1600/rings1600_msd_origins_full.dat'

# g1 or g3
g = 'g1'

# window parameters
windows = 40

# monomers
if('100' in fin): monomers = 100
if('200' in fin): monomers = 200
if('400' in fin): monomers = 400
if('800' in fin): monomers = 800
if('1600' in fin): monomers = 1600

# read data
t =  pypolymer.get_data(fin, 0)
g1 = pypolymer.get_data(fin, 1)
g3 = pypolymer.get_data(fin, 3)
total_points = len(t)

print 'first line = ', t[0], g1[0], g3[0]
print 'last line = ', t[-1], g1[-1], g3[-1]
print 'total points = ', total_points

bounds = np.logspace(2, 8, num=windows)
intervals = bounds.size - 1

x_ = t
y_ = g3 if(g == 'g3') else g1

# get correct starting index and trim list if necessary
start = 0
while (x_[start] < bounds[0]):
  start = start + 1           
  if(start >= total_points):
    sys.stdout.write('Bound error.\n')
    sys.exit(1)
if(start != 0):
  x_ = x_[start:]
  y_ = y_[start:]

jstart = 0
x = []
y = []
slopes = []
for i in range(intervals):
  lower_bound = bounds[i]
  upper_bound = bounds[i+1]
  t_seg = []
  p_seg = []
  for j in range(jstart, len(x_)):
    if(lower_bound <= x_[j] < upper_bound):                                                                                      
      t_seg.append(x_[j])                                                                                                         
      p_seg.append(y_[j])                                                                                                         
    else:
      if(len(t_seg) > 0 or x_[j] > upper_bound):
	jstart = j
	break
  if(len(t_seg) >= 2):
    xlog = np.log(np.array(t_seg))
    ylog = np.log(np.array(p_seg))
    m, b = np.polyfit(xlog, ylog, 1)
    a = math.log(t_seg[0],  10)
    b = math.log(t_seg[-1], 10)
    t_mid = 10**(0.5 * (a + b))
    slopes.append([t_mid, m])
    print '%d %d %d %.1f %.3f' % (t_seg[0], t_seg[-1], len(t_seg), t_mid, m)

outfile = 'slopes_N' + str(monomers) + '_' + g + '_W' + str(windows) + '.dat'
f = open(outfile, 'w')
for s in slopes:
  f.write('%.1f %.3f\n' % (s[0], s[1]))
f.close()

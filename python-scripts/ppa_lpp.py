#!/usr/bin/env python

"""This Python script computes the average contour length of
   a system of ring and linear polymers."""

import sys
import numpy as np
import glob
import pypolymer

# time step (tau)
dt = 0.01

blockfiles = glob.glob('ppa.*')
blockfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), blockfiles)
blockfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))

data_rings = []
data_linear = []
for file in blockfiles:

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

  sys.stdout.write('Working on ' + file + ' ...\n')

  re2_rings = {100:50.76, 200:88.75, 400:149.41, 800:242.41, 1600:401.73}
  re2_linear = {100:263.79, 200:538.94, 400:1095.33, 800:2164.92}

  if(file == blockfiles[0]):
    sys.stdout.write('(M_rings, M_linear, monomers) = (%d, %d, %d)\n' % (M_rings, M_linear, monomers))

  # compute the distance between two points in Euclidean space
  def rij(x1_, y1_, z1_, x2_, y2_, z2_):
    return ((x1_ - x2_)**2 + (y1_ - y2_)**2 + (z1_ - z2_)**2)**0.5

  # rings
  if(M_rings != 0):
    Lpp_rings = []
    for i in range(M_rings):
      clength = 0.0
      for j in range(monomers - 1):
        p1 = i * monomers + j
        p2 = p1 + 1
        clength = clength + rij(x[p1], y[p1], z[p1], x[p2], y[p2], z[p2])
      p1 = i * monomers + monomers - 1
      p2 = i * monomers
      clength = clength + rij(x[p1], y[p1], z[p1], x[p2], y[p2], z[p2])
      Lpp_rings.append(clength)

  # linear
  if(M_linear != 0):
    Lpp_linear = []
    for i in range(M_linear):
      clength = 0.0
      for j in range(monomers - 1):
        p1 = M_rings * monomers + i * monomers + j
        p2 = p1 + 1
        clength = clength + rij(x[p1], y[p1], z[p1], x[p2], y[p2], z[p2])
      Lpp_linear.append(clength)

  if(M_rings != 0):
    ave_Lpp_rings = np.mean(np.array(Lpp_rings))
    Ne_rings = 0.5 * monomers * re2_rings[monomers] / (0.5 * ave_Lpp_rings)**2
    data_rings.append((timestep * dt, ave_Lpp_rings, Ne_rings))

  if(M_linear != 0):
    ave_Lpp_linear = np.mean(np.array(Lpp_linear))
    Ne_linear = monomers * re2_linear[monomers] / ave_Lpp_linear**2
    data_linear.append((timestep * dt, ave_Lpp_linear, Ne_linear))

if(M_rings != 0):
  outfile = 'Lpp_' + str(M_rings) + '_' + str(M_linear) + '_N' + str(monomers) + '_rings.dat'
  f = open(outfile, 'w')
  f.write('# time (tau)  Lpp_rings  Ne_rings\n')
  for d in data_rings:
    f.write('%.3f %.3f %.3f\n' % (d[0], d[1], d[2]))
  f.close()
  sys.stdout.write(outfile + ' has been written.\n')

if(M_linear != 0):
  outfile = 'Lpp_' + str(M_rings) + '_' + str(M_linear) + '_N' + str(monomers) + '_linear.dat'
  f = open(outfile, 'w')
  f.write('# time (tau)  Lpp_linear  Ne_linear\n')
  for d in data_linear:
    f.write('%.3f %.3f %.3f\n' % (d[0], d[1], d[2]))
  f.close()
  sys.stdout.write(outfile + ' has been written.\n')

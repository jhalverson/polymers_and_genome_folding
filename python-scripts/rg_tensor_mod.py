"""This script computes the eigenvalue ratios of the
   gyration tensor and Rg squared of each chain for
   each configuration."""

# use sys for printing compatibility between Python 2 & 3
from __future__ import with_statement
import sys
import math
import os
import time
import numpy
import pypolymer

# case (a, b or c)
case = 'a'

# intermediate output ('yes' or 'no')
verbose = 'yes'

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('../../rings800a.*')

  # sort blockfiles list
  blockfiles.sort(lambda u, v: cmp(int(u[u.index('.', 5) + 1:]), int(v[v.index('.', 5) + 1:])))

  # work with subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M1600.0816000000']
 
else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

# loop over files
for file in blockfiles:

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

  r = []
  for i in range(total_particles):
    r.append((x[i], y[i], z[i]))

  if(file == blockfiles[0]):
    timeval_zero = timestep / 100
  timeval = timestep / 100

  # initialize list to store data
  L13_lst = []
  L23_lst = []
  RGsq_lst = []

  # loop over chains
  for i in xrange(M_rings):

    # initialize the tensor
    S = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    # compute center-of-mass
    rcm = [0.0, 0.0, 0.0]
    for j in xrange(monomers):
      p = i * monomers + j
      rcm[0] = rcm[0] + r[p][0]
      rcm[1] = rcm[1] + r[p][1]
      rcm[2] = rcm[2] + r[p][2]
    rcm[0] = rcm[0] / monomers
    rcm[1] = rcm[1] / monomers
    rcm[2] = rcm[2] / monomers

    # loop over components and then particles
    for m in [0, 1, 2]:
      for n in [0, 1, 2]:
        for j in xrange(monomers):
          p = i * monomers + j
          S[m][n] = S[m][n] + (r[p][m] - rcm[m]) * (r[p][n] - rcm[n])

    # normalize
    for m in [0, 1, 2]:
      for n in [0, 1, 2]:
	S[m][n] = S[m][n] / monomers

    # convert to a numpy array, compute eigenvalues and sort
    P = numpy.array(S)
    eigvals = numpy.linalg.eigvalsh(P)
    eigvals.sort()

    # compute ratios and radius of gyration squared
    L13 = eigvals[2] / eigvals[0]
    L23 = eigvals[1] / eigvals[0]
    RGsq = sum(eigvals)

    L13_lst.append(L13)
    L23_lst.append(L23)
    RGsq_lst.append(RGsq)

  # file L13
  s = ''
  for i in range(M_rings):
    s = s + '%12.4f ' % (L13_lst[i])
  fout = 'L13_all_' + str(monomers) + case + '.dat'
  with open(fout, 'a') as f:
    f.write('%10d %10d' % (timeval - timeval_zero, timeval) + ' ' + s + '\n')

  # file L23
  s = ''
  for i in range(M_rings):
    s = s + '%12.4f ' % (L23_lst[i])
  fout = 'L23_all_' + str(monomers) + case + '.dat'
  with open(fout, 'a') as f:
    f.write('%10d %10d' % (timeval - timeval_zero, timeval) + ' ' + s + '\n')

  # file RGsq
  s = ''
  for i in range(M_rings):
    s = s + '%12.4f ' % (RGsq_lst[i])
  fout = 'RGsq_all_' + str(monomers) + case + '.dat'
  with open(fout, 'a') as f:
    f.write('%10d %10d' % (timeval - timeval_zero, timeval) + ' ' + s + '\n')

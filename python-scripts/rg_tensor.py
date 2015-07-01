"""This script computes the time-averaged eigenvalues of the
   gyration tensor of each chain. It only looks at M_rings
   and ignores M_linear (which may be of interest)."""

import sys
import os
import glob
import numpy
import pypolymer

# case (a, b or c)
case = 'a'

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  blockfiles = glob.glob('../../rings1600a.*')
  blockfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), blockfiles)
  blockfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = blockfiles[::]

elif(blockfiles_method == 2):

  prefix = ''
  istr = 0
  iend = 10000000
  incr = 100000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append(prefix + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M1600.0495000000']
 
else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

# initialize eigenvalues
lambda1 = 0.0; lambda2 = 0.0; lambda3 = 0.0

# loop over files
for file in blockfiles:

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

  sys.stdout.write('Working on ' + file + ' ...\n')

  r = []
  for i in range(len(x)):
    r.append((x[i], y[i], z[i]))

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

    # following the convention of Vettorel et al. 2009
    lambda1 = lambda1 + eigvals[2]
    lambda2 = lambda2 + eigvals[1]
    lambda3 = lambda3 + eigvals[0]

# average eigenvalues
lambda1_ave = lambda1 / (M_rings * len(blockfiles))
lambda2_ave = lambda2 / (M_rings * len(blockfiles))
lambda3_ave = lambda3 / (M_rings * len(blockfiles))

# write data to standard output
sys.stdout.write('\n')
sys.stdout.write('lambda1 lambda2 lambda3: %.2f %.2f %.2f sigma * sigma\n' % (lambda1_ave, lambda2_ave, lambda3_ave))
sys.stdout.write('RG squared: %.2f sigma * sigma\n' % (lambda1_ave + lambda2_ave + lambda3_ave))
sys.stdout.write('asphericity: %4.2f sigma * sigma\n' % (lambda1_ave - 0.5 * (lambda2_ave + lambda3_ave)))
sys.stdout.write('lambda13: %.2f\n' % (lambda1_ave / lambda3_ave))
sys.stdout.write('lambda23: %.2f\n' % (lambda2_ave / lambda3_ave))

# write data to file
outfile = 'eigenvalues' + str(monomers) + case + '.dat'
fout = open(outfile, 'w')
fout.write('# lambda1 lambda2 lambda3 (sigma**2)\n')
fout.write('%8.2f %8.2f %8.2f\n\n' % (lambda1_ave, lambda2_ave, lambda3_ave))
fout.write('# RG squared (sigma**2)\n')
fout.write('%8.2f\n\n' % (lambda1_ave + lambda2_ave + lambda3_ave))
fout.write('# lambda13 lambda23\n')
fout.write('%8.2f %8.2f\n' % (lambda1_ave / lambda3_ave, lambda2_ave / lambda3_ave))
fout.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

"""This script computes the self-density of chains."""

import sys
import time
import os
import math
import glob
import pypolymer

case = ''

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 2

if(blockfiles_method == 1):

  prefix = '../../'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]
  if 0:
    maxfiles = 30
    blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))]

elif(blockfiles_method == 2):

  prefix = '../../'
  istr = 0
  iend = 10000000
  incr = 100000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append(prefix + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['../../']

else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit()

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                      rho_self                   \n')
sys.stdout.write('=================================================\n')
sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')
if(len(blockfiles) > 2):
  sys.stdout.write('  blockfiles[0] = ' + blockfiles[0] + '\n')
  sys.stdout.write('  blockfiles[1] = ' + blockfiles[1] + '\n')
  sys.stdout.write('  blockfiles[2] = ' + blockfiles[2] + '\n')
  sys.stdout.write('  ...' + '\n')
  sys.stdout.write('  blockfiles[-3] = ' + blockfiles[-3] + '\n')
  sys.stdout.write('  blockfiles[-2] = ' + blockfiles[-2] + '\n')
  sys.stdout.write('  blockfiles[-1] = ' + blockfiles[-1] + '\n\n')

# initialize the histogram
dr = 0.1
max_r = 20.0
max_bin = int(max_r / dr)
hist = []
for i in range(max_bin):
  hist.append(0)

for file in blockfiles:

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

  if(file == blockfiles[0]):
    total_particles1st = total_particles
    sys.stdout.write('System parameters:\n')
    sys.stdout.write('  M_rings = ' + str(M_rings) + '\n')
    sys.stdout.write('  M_linear = ' + str(M_linear) + '\n')
    sys.stdout.write('  N = ' + str(monomers) + '\n')
    sys.stdout.write('  N/2 = ' + str(monomers_half) + '\n\n')
  if(total_particles1st != total_particles):
    sys.stderr.write('Change in total_particles (' + str(total_particles) + ').\n')
    sys.exit()

  sys.stdout.write('Working on ' + file + '\n')

  for i in range(M_rings + M_linear):
    xcm = 0.0
    ycm = 0.0
    zcm = 0.0
    for j in range(monomers):
      p = i * monomers + j
      xcm = xcm + x[p]
      ycm = ycm + y[p]
      zcm = zcm + z[p]
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers

    for j in range(monomers):
      p = i * monomers + j
      rij = ((xcm - x[p])**2 + (ycm - y[p])**2 + (zcm - z[p])**2)**0.5
      ibin = int(rij / dr)
      if(ibin >= max_bin): print 'ERROR: max_bin too small'
      hist[ibin] = hist[ibin] + 1

# normalize
outfile = 'rho' + case + '.dat'
fout = open(outfile, 'w')
for i in range(max_bin):
  inner = i * dr
  outer = inner + dr
  mid = 0.5 * (inner + outer)
  vol = (4.0 * math.pi / 3.0) * (outer**3 - inner**3)
  rho = hist[i] / (vol * len(blockfiles) * (M_linear + M_rings))
  fout.write('%.4f %.4f %.4f %.4f\n' % (inner, outer, mid, rho))
fout.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

"""This Python script computes the mean-square internal distances
   for a system of ring polymers."""

import sys
import math
import glob
import pypolymer
import numpy as np

# show result in Matplotlib
mpl = False

# generation of file names: (1) glob, (2) loop, (3) manual or (4) log scale
blockfiles_method = 2

if(blockfiles_method == 1):

  prefix = 'rings400a.'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]
  if 0:
    maxfiles = 30
    blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))] = blockfiles[::]

elif(blockfiles_method == 2):

  istr = 2000000
  iend = 2000000
  incr = 2000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('../sqt/stripped/' + str(i) + '.pos')

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M400.end']

elif(blockfiles_method == 4):

  # prefix
  prefix = '../../rings200_26.'

  # number of files per decade
  n = 5

  # exponent of first and last decade (time steps)
  start_exp = 5
  end_exp = 10
 
  # base
  b = math.e

  number_of_decades = end_exp - start_exp
  start = math.log(10**start_exp)
  end = math.log(10**end_exp)
  step = (end - start) / (n * number_of_decades - 1)
  uniform_time = []
  for i in range(n * number_of_decades):
    uniform_time.append(b**(start + step * i))
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.index('.', 5) + 1:].isdigit(), existingfiles)
  existingfiles = [int(u[u.index('.', 5) + 1:]) for u in existingfiles]
  existingfiles.sort()
  uniform_time = [pypolymer.find_closest_existing_file(existingfiles, u) for u in uniform_time]
  uniform_time = pypolymer.unique(uniform_time)
  blockfiles = [prefix + str(u) for u in uniform_time]

else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit()

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                      [d(s)]^2                   \n')
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
    dsq = np.zeros(monomers_half + 1, dtype=np.float64)
    dsq_dsqave = np.zeros(monomers_half + 1, dtype=np.float64)
    ct = np.zeros(monomers_half + 1, dtype=np.int32)
  if(total_particles1st != total_particles):
    sys.stderr.write('Change in total_particles (' + str(total_particles) + ').\n')
    sys.exit()

  sys.stdout.write('Working on ' + file + '\n')

  for i in range(M_rings):
    for j in range(monomers - 1):
      p1 = i * monomers + j
      xp1 = x[p1]
      yp1 = y[p1]
      zp1 = z[p1]
      for k in range(j + 1, monomers):
        p2 = i * monomers + k
        rsq = (xp1 - x[p2])**2 + (yp1 - y[p2])**2 + (zp1 - z[p2])**2
        s = k - j
        if(s > monomers_half): s = monomers - s
        dsq[s] = dsq[s] + rsq
        ct[s] = ct[s] + 1

dsq = dsq / ct
for s in range(1, monomers_half + 1):
  dsq_dsqave[s] = dsq[s] * (1.0 / s + 1.0 / (monomers - s))

outfile = 'dsq' + str(monomers) + '.dat'
f = open(outfile, 'w')
f.write('# s  dsq(s)/sigma^2  dsq(s)/<dsq>\n')
for s in range(1, monomers_half + 1):
  f.write('%d %.3f %.3f\n' % (s, dsq[s], dsq_dsqave[s]))
f.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

if mpl:
  import matplotlib.pyplot as plt
  plt.semilogx(dsqdsq_ave, 'b-')
  plt.xlim(1, 1000)
  plt.ylim(0.1, 3)
  plt.show()

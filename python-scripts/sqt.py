"""This Python script computes the single-chain dynamic structure
   factor of a polymer melt."""

import sys
import os
import time
import math
import glob
import pypolymer
import numpy as np
import cmath
import random

# integration time step (tau)
delta_t = 0.01

# list of q magnitudes (sigma^-1)
q_magn = np.array([0.2, 0.4, 0.6, 0.8, 1.0])

# number of random vectors
random_vectors_per_q_magn = 10

# job case
case = ''

# displacements for time origins (range() may be used)
origins = []
origins = map(int, origins)

# seed the RNG
random.seed(10)

# generation of file names: (1) glob, (2) loop, (3) manual or (4) log scale
blockfiles_method = 4

if(blockfiles_method == 1):

  prefix = 'rings400a.'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]
  if 0:
    maxfiles = 30
    blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))]

elif(blockfiles_method == 2):

  prefix = '../../rings200_26.'
  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append(prefix + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['rings400a.1000000000']

elif(blockfiles_method == 4):

  # prefix
  prefix = '../../rings200_26.'

  # maximum number of files (script may use less)
  n = 100

  # exponent of first and last decade (time steps)
  start_exp = 5
  end_exp = 10
 
  step_exp = float(end_exp - start_exp) / (n - 1)
  uniform_time = [10**(start_exp + step_exp * u) for u in range(n)]
  existingfiles = glob.glob(prefix + '*')
  if(existingfiles == []): sys.stderr.write('ERROR: No files found.\n'); sys.exit(1)
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles = [int(u[u.rindex('.') + 1:]) for u in existingfiles]
  existingfiles.sort()
  uniform_time = [pypolymer.find_closest_existing_file(existingfiles, u) for u in uniform_time]
  uniform_time = pypolymer.unique(uniform_time)
  blockfiles = [prefix + str(u) for u in uniform_time]

else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit()

if(len(blockfiles) == 0):
  sys.stderr.write('ERROR: blockfiles is empty\n')
  sys.exit()

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                       S(q,t)                    \n')
sys.stdout.write('=================================================\n')
sys.stdout.write('\n\nTotal number of files in first set: ' + str(len(blockfiles)) + '\n')
if(len(blockfiles) > 2):
  sys.stdout.write('  blockfiles[0] = ' + blockfiles[0] + '\n')
  sys.stdout.write('  blockfiles[1] = ' + blockfiles[1] + '\n')
  sys.stdout.write('  blockfiles[2] = ' + blockfiles[2] + '\n')
  sys.stdout.write('  ...' + '\n')
  sys.stdout.write('  blockfiles[-3] = ' + blockfiles[-3] + '\n')
  sys.stdout.write('  blockfiles[-2] = ' + blockfiles[-2] + '\n')
  sys.stdout.write('  blockfiles[-1] = ' + blockfiles[-1] + '\n\n')

def getRand(b_):
  """Returns a vector of magnitude b_ from a uniform distribution on the surface of a sphere."""
  x1 = x2 = 1.0
  while(x1**2 + x2**2 >= 1.0):
    x1 = 2.0 * random.random() - 1.0
    x2 = 2.0 * random.random() - 1.0
  x = 2.0 * x1 * (1 - x1**2 - x2**2)**0.5
  y = 2.0 * x2 * (1 - x1**2 - x2**2)**0.5
  z = 1.0 - 2.0 * (x1**2 + x2**2)
  return (x * b_, y * b_, z * b_)

# create list of file sets
fset = [blockfiles]

if(len(origins) != 0):
  # store existing files
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)

  # build file sets for the different time origins
  for i in range(len(origins)):
    fs = []
    disp = origins[i]
    for file in blockfiles:
      suffix = int(file[file.rindex('.') + 1:])
      fname = prefix + str(suffix + disp)
      fs.append(fname)
    fs = list(set(fs) & set(existingfiles))
    fs.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
    fset.append(fs)

sys.stdout.write('Number of file sets = ' + str(len(fset)) + '\n')
for i in range(len(fset)):
  sys.stdout.write('  number of files in set ' + str(i) + ' = ' + str(len(fset[i])) + '\n')

# initialize dictionaries
S_ave = []
for i in range(q_magn.size):
  S_ave.append({})
cntr = {}

# for each file set
for files in fset:

  xi = []; yi = []; zi = []
  sys.stdout.write('\n\nStarting new file set ...\n')

  qvec_lst = []
  for magn in q_magn:
    tmp = []
    for i in range(random_vectors_per_q_magn):
      tmp.append(getRand(magn))
    qvec_lst.append(tmp)

  # for each file in the file set
  for file in files:

    fmanager = pypolymer.read_blockfile(file)
    timestep = fmanager[0]
    total_particles = fmanager[1]
    M_rings = fmanager[2]; M_linear = fmanager[3]
    monomers = fmanager[4]; monomers_half = fmanager[5]
    x = fmanager[6]; y = fmanager[7]; z = fmanager[8]
    M_rings = 113
    M_linear = 113

    if(files == fset[0] and file == files[0]):
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

    # use first file of each file set to get initial values
    if(file == files[0]):
      # store reference time
      time_zero = timestep
      if(timestep != 0):
	sys.stdout.write('WARNING: time is not 0 for initial file (' + str(timestep) +  ').\n\n')
      # store monomer positions at t=0
      for i in range(M_rings):
	for j in range(monomers):
	  p = i * monomers + j
	  xi.append(x[p])
	  yi.append(y[p])
	  zi.append(z[p])
      # precompute t=0 sums
      isum = []
      for i in range(M_rings):
        tmp1 = []
        for qvec in qvec_lst:
          tmp2 = []
          for q in qvec:
	    sum = complex(0.0, 0.0)
	    for m in range(monomers):
              p = i * monomers + m
              sum = sum + cmath.exp(-complex(0, 1) * (q[0] * xi[p] + q[1] * yi[p] + q[2] * zi[p]))
	    tmp2.append(sum)
          tmp1.append(tmp2)
        isum.append(tmp1)

    # loop over rings
    S = np.zeros(q_magn.size)
    for i in range(M_rings):
      for j, qvec in enumerate(qvec_lst):
        for k, q in enumerate(qvec):
          sum = complex(0.0, 0.0)
          for m in range(monomers):
            p = i * monomers + m
            sum = sum + cmath.exp(complex(0, 1) * (q[0] * x[p] + q[1] * y[p] + q[2] * z[p]))
          prod = isum[i][j][k] * sum
          S[j] = S[j] + prod.real

    # store time difference
    dt = timestep - time_zero

    if(not cntr.has_key(dt)):
      cntr[dt] = 0
      for i in range(q_magn.size):
        S_ave[i][dt] = 0.0
 
    cntr[dt] = cntr[dt] + 1
    for i in range(q_magn.size):
      S_ave[i][dt] = S_ave[i][dt] + S[i] / (M_rings * monomers * random_vectors_per_q_magn)

# normalize data
data = []
for dt in sorted(cntr.keys()):
  s = ''
  for i in range(q_magn.size):
    S_ave[i][dt] = S_ave[i][dt] / cntr[dt]
    s = s + '%.3f  ' % (S_ave[i][dt])
  data.append((dt, s))

# write data to file
outfile = 'Sqt' + str(monomers) + case + '.dat'
fout = open(outfile, 'w')
fout.write('# ' + time.asctime() + '\n')
fout.write('# ' + os.getcwd() + '\n')
fout.write('# M_rings  M_linear  monomers\n')
fout.write('# ' + str(M_rings) + '  ' + str(M_linear) + '  ' + str(monomers) + '\n')
fout.write('# S(q,t)\n')
s = '# time (tau)  q='
for i in range(q_magn.size):
  s = s + '%.1f  ' % (q_magn[i])
fout.write('%s(sigma^-1)\n' % s)
for item in data:
  fout.write('%9d %s\n' % (delta_t * item[0], item[1]))
fout.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

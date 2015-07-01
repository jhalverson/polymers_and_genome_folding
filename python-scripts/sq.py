"""This Python script computes the single-chain static structure factor
   for a system of ring polymers."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import glob
import pypolymer
import numpy as np
import random

# magnitude of q vectors (10^-2, 10^1)
q_magn = np.logspace(-2, 1, num=80)

# numpy array to store S(q)
S = np.zeros(q_magn.size)

# number of random vectors for each |q|
random_vectors_per_q_magn = 5

# show result in Matplotlib
mpl = True

# generation of file names: (1) glob, (2) loop, (3) manual or (4) log scale
blockfiles_method = 3

if(blockfiles_method == 1):

  prefix = 'rings400a.'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.index('.', 5) + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.index('.', 5) + 1:]), int(v[v.index('.', 5) + 1:])))
  blockfiles = existingfiles[::]
  if 0:
    maxfiles = 30
    blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))] = blockfiles[::]

elif(blockfiles_method == 2):

  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('../../rings200_26.' + str(i))

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
sys.stdout.write('                        S(q)                     \n')
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

def getRand(b_):
  """Returns uniformly distributed vectors on the surface of a sphere with radius b_."""
  x1 = x2 = 1.0
  while(x1**2 + x2**2 >= 1.0):
    x1 = 2.0 * random.random() - 1.0
    x2 = 2.0 * random.random() - 1.0
  x = 2.0 * x1 * (1 - x1**2 - x2**2)**0.5
  y = 2.0 * x2 * (1 - x1**2 - x2**2)**0.5
  z = 1.0 - 2.0 * (x1**2 + x2**2)
  return (x * b_, y * b_, z * b_)

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

  for i in range(M_rings):
    for j, magn in enumerate(q_magn):
      for k in range(random_vectors_per_q_magn):
	q = getRand(magn)
	sin_sum = 0.0
	cos_sum = 0.0
	for l in range(monomers):
	  p = i * monomers + l
	  q_dot_r = q[0] * x[p] + q[1] * y[p] + q[2] * z[p]
	  sin_sum = sin_sum + math.sin(q_dot_r)
	  cos_sum = cos_sum + math.cos(q_dot_r)
	S[j] = S[j] + sin_sum**2 + cos_sum**2

#   alternative version (must import cmath)
#   for i in range(M_rings):
#    for j, magn in enumerate(q_magn):
#      for k in range(random_vectors_per_q_magn):
#        q = getRand(magn)
#        sum = complex(0.0, 0.0)
#        for l in range(monomers):
#          p = i * monomers + l
#          q_dot_r = q[0] * x[p] + q[1] * y[p] + q[2] * z[p]
#          sum = sum + cmath.exp(complex(0, 1) * q_dot_r)
#        S[j] = S[j] + sum * sum.conjugate()

S = S / (M_rings * random_vectors_per_q_magn * monomers * len(blockfiles))

outfile = 'Sq' + str(monomers) + '.dat'
f = open(outfile, 'w')
f.write('# q (sigma^-1)  S(q)\n')
for i in range(q_magn.size):
  f.write('%.3f %.3f\n' % (q_magn[i], S[i]))
f.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

if mpl:
  import matplotlib.pyplot as plt
  plt.loglog(q_magn, S, 'bo')
  plt.loglog(q_magn, S, 'k-', lw=0.5)
  plt.xlim(0.01, 1e1)
  plt.ylim(0.1, 2e3)
  plt.show()

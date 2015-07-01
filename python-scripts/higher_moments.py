"""This script computes Rg^2, Rg^4, Rg^6, Rg^8 and Rg^10 for a system
   of rings and/or linear chains. The average values and
   standard deviations are reported. Coordinates are assumed
   to be unwrapped."""

# use sys for printing compatibility between Python 2 & 3
import sys
import time
import os
import math
import glob
import pypolymer
import numpy as np

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
    blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))]

elif(blockfiles_method == 2):

  istr = 1500000000
  iend = 2000000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('../../200/' + str(i) + '.pos')

elif(blockfiles_method == 3):

  blockfiles = ['rings200_26.0']

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
sys.stdout.write('                      Rg2 Re2                    \n')
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

# initialize the lists (time, Rgsq, Resq)
rg2 = []
rg4 = []
rg6 = []
rg8 = []
rg1 = []

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

    rad2 = 0.0
    rad4 = 0.0
    rad6 = 0.0
    rad8 = 0.0
    rad1 = 0.0
    for j in range(monomers):
      p = i * monomers + j
      tmp = xcm**2 + ycm**2 + zcm**2 - 2.0 * (x[p] * xcm + y[p] * ycm + z[p] * zcm) + x[p]**2 + y[p]**2 + z[p]**2
      rad2 = rad2 + tmp
      rad4 = rad4 + tmp**2
      rad6 = rad6 + tmp**3
      rad8 = rad8 + tmp**4
      rad1 = rad1 + tmp**5
    rad2 = rad2 / monomers
    rad4 = rad4 / monomers
    rad6 = rad6 / monomers
    rad8 = rad8 / monomers
    rad1 = rad1 / monomers
    rg2.append(rad2)
    rg4.append(rad4)
    rg6.append(rad6)
    rg8.append(rad8)
    rg1.append(rad1)

x2 = np.array(rg2)
print '<Rg^2>/sigma^2:', np.mean(x2)
x4 = np.array(rg4)
print '<Rg^4>/sigma^4:', np.mean(x4)
x6 = np.array(rg6)
print '<Rg^6>/sigma^6:', np.mean(x6)
x8 = np.array(rg8)
print '<Rg^8>/sigma^6:', np.mean(x8)
x1 = np.array(rg1)
print '<Rg^10>/sigma^10:', np.mean(x1)

print '<Rg^(4)>/Rg^(2):', np.mean(x4)**(1/4.0) / np.mean(x2)**(1/2.0)
print '<Rg^(6)>/Rg^(2):', np.mean(x6)**(1/6.0) / np.mean(x2)**(1/2.0)
print '<Rg^(8)>/Rg^(2):', np.mean(x8)**(1/8.0) / np.mean(x2)**(1/2.0)
print '<Rg^(10)>/Rg^(2):', np.mean(x1)**(1/10.0) / np.mean(x2)**(1/2.0)


print '<Rg^4>/<Rg^2>^2 = ', np.mean(x4) / np.mean(x2)**2
print '<Rg^6>/<Rg^2>^3 = ', np.mean(x6) / np.mean(x2)**3
print '<Rg^8>/<Rg^2>^4 = ', np.mean(x8) / np.mean(x2)**4
print '<Rg^10>/<Rg^2>^5 = ', np.mean(x1) / np.mean(x2)**5
print '<Rg^8>/<Rg^4>^2 = ', np.mean(x8) / np.mean(x4)**2

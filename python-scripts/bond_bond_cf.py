"""This Python script computes the bond-bond correlation function
   for linear chains."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import pypolymer
import numpy as np

# number of chains M
chains = 400

# monomers
monomers = 800

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('*.pos')
  #blockfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), blockfiles)
  #blockfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))

  # work with a subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 1000000
  iend = 2000000000
  incr = 1000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('../sqt/stripped/' + str(i) + '.pos')

elif(blockfiles_method == 3):

  blockfiles = ['rings400a.5000000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit()

sys.stdout.write('Total number of files: ' + str(len(blockfiles)) + '\n\n')

good_chains = 0
total = np.zeros(monomers - 1, dtype=np.float64)

# loop over block files
for file in blockfiles:

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]




  monomers = 800




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

  # begin calculation
  for i in range(chains):

    chain_good = True
    for j in range(monomers - 1):
      p1 = i * monomers + j
      p2 = p1 + 1
      dsq = (x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1]- z[p2])**2
      rij = dsq**0.5
      if(rij < 0.8 or rij > 1.2):
        print 'PROBLEM:', i, j, j + 1, rij
        chain_good = False
  
    if(chain_good):

      good_chains = good_chains + 1

      r21 = []
      for j in range(monomers - 1):
        p1 = i * monomers + j
        p2 = p1 + 1
        r21.append([x[p2] - x[p1], y[p2] - y[p1], z[p2] - z[p1]])

      t_run = len(r21)
      cf = []
      for t in range(monomers - 1):
	t_max = t_run - t
	sum_aa = 0.0
	for t0 in range(t_max):
	  sum_aa = sum_aa + r21[t0][0] * r21[t0 + t][0] \
			  + r21[t0][1] * r21[t0 + t][1] \
			  + r21[t0][2] * r21[t0 + t][2]
	cf.append(sum_aa / t_max)
      #cf.insert(0, 0)
      total = total + np.array(cf)
total = total / good_chains

print good_chains, chains * len(blockfiles)

f = open('ps1.dat', 'w')
for i in range(monomers - 1):
  f.write('%d %.6e\n' % (i+1, total[i] / total[0]))
f.close()

sys.stdout.write('Output has been written to disk.\n')

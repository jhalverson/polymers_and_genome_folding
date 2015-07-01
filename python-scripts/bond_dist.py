"""This Python script computes the distribution of bond lengths
   for a polymer melt of linear or rings."""

import sys
import os
import time
import numpy as np
import pypolymer

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 2
mpl = True

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('linear_npt.*')
  blockfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), blockfiles)
  blockfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))

  # work with a subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 500000
  iend = 1000000
  incr = 10000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('linear_npt.' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['rings400a.5000000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit()

sys.stdout.write('Total number of files: ' + str(len(blockfiles)) + '\n\n')

bonds = []
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
    if(len(blockfiles) * (M_rings + M_linear) > 100000):
      print 'This will produce more than 1e5 values. Exiting ...'
      sys.exit()
  if(total_particles1st != total_particles):
    sys.stderr.write('Change in total_particles (' + str(total_particles) + ').\n')
    sys.exit()

  sys.stdout.write('Working on ' + file + '\n')

  # rings
  for i in range(M_rings):
    for j in range(monomers - 1):
      p1 = i * monomers + j
      p2 = p1 + 1
      rij = ((x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1] - z[p2])**2)**0.5
      bonds.append(rij)
    p1 = i * monomers
    p2 = p1 + monomers - 1
    rij = ((x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1] - z[p2])**2)**0.5
    bonds.append(rij)

  # linear
  for i in range(M_linear):
    for j in range(monomers - 1):
      p1 = M_rings * monomers + i * monomers + j
      p2 = p1 + 1
      rij = ((x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1] - z[p2])**2)**0.5
      bonds.append(rij)

# write data to file
outfile = 'bonds_Mrings' + str(M_rings) + '_Mlinear' + str(M_linear) + '_N' + str(monomers) + '.dat'
fout = open(outfile, 'w')
fout.write('# ' + time.asctime() + '\n')
fout.write('# ' + os.getcwd() + '\n')
fout.write('# M_rings M_linear monomers\n')
fout.write('# ' + str(M_rings) + ' ' + str(M_linear) + ' ' + str(monomers) + '\n')
for i in range(len(bonds)):
  fout.write('%6.4f\n' % bonds[i])
fout.close()

if(mpl):

  import matplotlib.pyplot as plt

  plt.rcParams['text.usetex'] = True
  plt.rcParams['legend.numpoints'] = 1
  plt.rcParams['axes.linewidth'] = 0.5
  plt.rcParams['axes.titlesize'] = 12
  plt.rcParams['axes.labelsize'] = 10
  plt.rcParams['xtick.labelsize'] = 8
  plt.rcParams['ytick.labelsize'] = 8
  plt.rcParams['font.size'] = 8

  w = 8.0 / 2.54
  h = (3.0 / 4.0) * w

  bonds_ave = np.mean(bonds)
  bonds_std = np.std(bonds)
  ave_str = '%.3f' % bonds_ave
  std_str = '%.3f' % bonds_std
  print '<b>=', ave_str, ' std=', std_str

  fig = plt.figure(1, figsize=(w, h))
  ax1 = fig.add_subplot(1, 1, 1)
  plt.hist(bonds, bins=50, normed=True)
  plt.xlabel(r'$b/\sigma$', fontsize=10)
  plt.ylabel(r'$P(b)$', fontsize=10)
  #plt.title(r'$\mathrm{Rings}$ $N = 800$', fontsize=10)
  plt.text(1.01, 11, r'$\langle b \rangle =' + ave_str + '\pm' + std_str + '~\sigma$', fontsize=8)

  upper_right = (0.85, 0.85)
  from_left = 0.20
  lower_left = (from_left, 1.0 - ((upper_right[0] - from_left) + (1.0 - upper_right[1])))

  fig.subplots_adjust(left=lower_left[0], right=upper_right[0], bottom=lower_left[1], \
                      top=upper_right[1], wspace=0.0, hspace=0.0)

  if 1:
    for ax in fig.axes:
      plt.setp(ax.patch, color=(0.9, 0.9, 0.9))

  plt.savefig('bond_dist.png', dpi=200)
  os.system('xv bond_dist.png&')

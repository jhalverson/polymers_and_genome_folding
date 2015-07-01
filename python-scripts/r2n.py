"""This Python script computes the mean-square internal distances
   for the linear chains of a ring/linear blend."""

# use sys for printing compatibility between Python 2 & 3
import sys
import time
import os
import math
import glob
import pypolymer

mpl = True
if mpl:
  import matplotlib.pyplot as plt

# generation of file names: (1) glob, (2) loop, (3) manual or (4) log scale
blockfiles_method = 2

if(blockfiles_method == 1):

  prefix = '../configs/rings200_6.'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]
  if 0:
    maxfiles = 30
    blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))]

elif(blockfiles_method == 2):

  istr = 3000000000
  iend = 3800000000
  incr = 500000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('../configs/rings200_6.' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['rings400_6_slow_insert.2450000']

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
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles = [int(u[u.rindex('.') + 1:]) for u in existingfiles]
  existingfiles.sort()
  uniform_time = [pypolymer.find_closest_existing_file(existingfiles, u) for u in uniform_time]
  uniform_time = pypolymer.unique(uniform_time)
  blockfiles = [prefix + str(u) for u in uniform_time]

else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit()

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                   r2(n)/n vs. n                 \n')
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

# get value of monomers
file = blockfiles[0]
fmanager = pypolymer.read_blockfile(file)
monomers = fmanager[4]

# initialization
norm = []
aveDistsq = []
for m in range(monomers):
  norm.append(0)
  aveDistsq.append(0.0)

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

  for k in range(M_linear):
    kstart = M_rings * monomers + k * monomers
    kend = kstart + monomers
    for i in range(kstart, kend - 1):
      xi = x[i]
      yi = y[i]
      zi = z[i]
      for j in range(i + 1, kend):
        m = j - i
        dsq = (xi - x[j])**2 + (yi - y[j])**2 + (zi - z[j])**2
        aveDistsq[m] = aveDistsq[m] + dsq / m
        norm[m] = norm[m] + 1

for m in range(1, monomers):
  aveDistsq[m] = aveDistsq[m] / norm[m]

print 'end-to-end distance squared', (monomers - 1) * aveDistsq[monomers - 1]

f = open('msid_' + str(M_rings) + '_' + str(M_linear) + '_N' + str(monomers) + '.dat', 'w')
for m in range(1, monomers):
  f.write('%d %10.5f\n' % (m, aveDistsq[m]))
f.close()

if(mpl):

  plt.rcParams['xtick.labelsize'] = 8
  plt.rcParams['ytick.labelsize'] = 8
  plt.rcParams['legend.numpoints'] = 1
  plt.rcParams['text.usetex'] = True
  plt.rcParams['axes.linewidth'] = 0.5

  fin = '/people/thnfs/homes/halvers/research/contaminants/r2ij_limiting_distribution/dist_best_combined.dat'
  x = pypolymer.get_data(fin, 0)
  y = pypolymer.get_data(fin, 1)

  w = 8.0 / 2.54
  h = (3.0 / 4.0) * w

  fig = plt.figure(1, figsize=(w, h))
  plt.semilogx(x, y, 'r-', label='ideal')
  plt.semilogx(range(0, monomers), aveDistsq, 'b-')
  plt.ylim(0, 3.5)
  plt.xlabel(r'$s$', fontsize=10)
  plt.ylabel(r'$[d(s)]^2/s$', fontsize=10)
  plt.title(r'$M_{\mathrm{rings}} = 200$, $M_{\mathrm{linear}} = 6$, $N = 200$, $\phi_{\mathrm{linear}} = 0.03$', ha='center', fontsize=8)
  lg = plt.legend(loc='lower left', prop={'size':8}, borderaxespad=1, handletextpad=0.5)
  lg.get_frame().set_linewidth(0.25)

  upper_right = (0.85, 0.85)
  from_left = 0.17
  lower_left = (from_left, 1.0 - ((upper_right[0] - from_left) + (1.0 - upper_right[1])))
  fig.subplots_adjust(left=lower_left[0], right=upper_right[0], bottom=lower_left[1], top=upper_right[1], wspace=0.0, hspace=0.0)

  if 1:
    for ax in fig.axes:
      plt.setp(ax.patch, color=(0.9, 0.9, 0.9))

  outfile = 'msid.png'
  plt.savefig(outfile, dpi=200)
  os.system('xv ' + outfile + '&\n')

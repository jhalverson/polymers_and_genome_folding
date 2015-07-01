"""This script computes the radius of gyration of each
   chain as a function of time."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os
import pypolymer

# all systems have two hundred chains
chains = 200

# intermediate output ('yes' or 'no')
verbose = 'yes'

# load data into Gnuplot ('yes' or 'no')
gnuplot_load = 'no'

# generation of file names: (1) glob, (2) loop, (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('*000000')
  blockfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), blockfiles)
  blockfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))

elif(blockfiles_method == 2):

  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RG_N200M1600.0361200000', 'RG_N200M1600.0361300000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

# which chains should be written to file: (1) all, (2) every nth, (3) only a few
chains2file = 1

if(chains2file == 1):
  # all chain data will be written
  indices = range(chains)
elif(chains2file == 2):
  # every nth chain
  nth = 20
  indices = range(0, chains, nth)
elif(chains2file == 3):
  # only these chains will have their data written
  indices = [0, 37, 40]
else:
  sys.stdout.write('Value of chains2file is not valid.\n\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

# initialize the list
time_rg = []

# loop over files
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

  # initialize the sublist
  time_rg_sub = []

  sys.stdout.write('Working on '+ file + '\n')

  # store the time value (tau)
  time_rg_sub.append(timestep)

  # loop over chains
  for i in range(chains):

    # compute center-of-mass
    xcm = 0.0; ycm = 0.0; zcm = 0.0
    for j in range(monomers):
      p = i * monomers + j
      xcm = xcm + x[p]
      ycm = ycm + y[p]
      zcm = zcm + z[p]
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers

    # sum the squares
    radsq = 0.0
    for j in range(monomers):
      p = i * monomers + j
      radsq = radsq + (x[p] - xcm)**2 + (y[p] - ycm)**2 + (z[p] - zcm)**2

    RG = math.sqrt(radsq / monomers)
    time_rg_sub.append(RG)

  # append sublist to list
  time_rg.append(time_rg_sub)

# write out the data
outfile = 'rg_history' + str(monomers) + '.dat'
fout = open(outfile, 'w')

# create comment line
s = ''
for index in indices:
  s = s + 'i=%-3d ' % (index)

# write comment line
fout.write('   #tau    '  + s + ' \n')

# write out data
for i in range(len(blockfiles)):
  s = ''
  for index in indices:
    s = s + '%5.2f ' % (time_rg[i][index + 1])
  time_fmt = '%10d ' % (time_rg[i][0])
  fout.write(time_fmt + s + '\n')
fout.close()
sys.stdout.write(outfile + ' has been written to disk (RG vs. time data file).\n')

# write a Gnuplot input file
gnufile = 'rg_history' + str(monomers) + '.plt'
f = open(gnufile, 'w')
f.write('set xlabel \"TIME (TAU)\" \n')
f.write('set ylabel \"RG (SIGMA)\" \n')
f.write("""plot \\\n""")
for i, j in enumerate(indices[0:len(indices) - 1]):
  f.write('\"' + outfile + '\" using 1:' + str(i + 2) + ' with lines title \"' + str(j) + '\", \\\n')
i = len(indices)
j = indices[-1]
f.write('\"' + outfile + '\" using 1:' + str(i + 1) + ' with lines title \"' + str(j) + '\" \n')
f.close()
sys.stdout.write(gnufile + ' has been written to disk (Gnuplot input file).\n')

if(gnuplot_load == 'yes'):
  cmd = 'gnuplot -persist ' + gnufile
  os.system(cmd)

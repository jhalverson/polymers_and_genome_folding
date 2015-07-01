"""This script computes the corrected g1, g2 and g3 by
   removing a constant center-of-mass velocity from the
   system.
"""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os

# all systems have two hundred chains
chains = 200

# define drift velocity
v = (1.24734e-05, 8.89521e-06, -8.27691e-06)

# intermediate output ('yes' or 'no')
verbose = 'yes'

# load data into Gnuplot ('yes' or 'no')
gnuplot_load = 'no'

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('*0000')

  # sort blockfiles list
  blockfiles.sort(lambda u, v: cmp(int(u[u.index('.') + 1:]), int(v[v.index('.') + 1:])))

  # work with a subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 10000000
  iend = 22500000
  incr = 500000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M800a_dpd.00000000', 'RGEQ_N200M800a_dpd.0100000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

# initialize list to store the corrected data
time_g123 = []

#initialize lists to store initial coordinates and COM positions
xi = []; yi = []; zi = []
com_zero = []

# loop over files
for file in blockfiles:

  # open file
  f = open(file)

  # find line where total particles is found
  line = ''
  while not '{n_part ' in line:
    line = f.readline()
  sline = line.replace('}','')
  total_particles = int(sline.split()[1])
  monomers = total_particles / chains

  # safety checks
  if(file == blockfiles[0]):
    monomers1st = monomers
    sys.stdout.write('chains = ' + str(chains) + ', monomers = ' + str(monomers) + '\n\n')
  if(monomers != monomers1st):
    sys.stdout.write('ERROR: monomers is not equal to monomers1st\n')
    sys.exit(1)

  sys.stdout.write('Working on '+ file + '\n')

  # find line where time is found
  line = ''
  while not '{time ' in line:
    line = f.readline()
  sline = line.replace('}','')
  timeval = int(eval(sline.split()[1]))

  # check that time matches file extension (tau is 100 dt)
  exten = int(file[file.find('.') + 1:]) / 100
  if(exten != timeval):
    sys.stdout.write('WARNING: Mismatch between time and file name ...\n')
    sys.stdout.write('WARNING: Using extension value of ' + str(exten) + \
                     ' instead of ' + str(timeval) + '\n\n') 
    timeval = exten

  # find line where coordinates begin
  line = ''
  while not '{particles ' in line:
    line = f.readline()

  # store the corrected coordinates
  x = []; y = []; z = []
  for i in range(total_particles):
    sline = f.readline().split()
    x.append(float(sline[1]) - v[0] * timeval)
    y.append(float(sline[2]) - v[1] * timeval)
    z.append(float(sline[3]) - v[2] * timeval)

  # use first file to get initial values
  if(file == blockfiles[0]):
    if(timeval != 0):
      sys.stdout.write('ERROR: time is not 0 for initial file.\n\n')
      sys.exit(1)
    # loop over chains
    for i in range(chains):
      # compute center-of-mass and store coordinates
      xicm = 0.0; yicm = 0.0; zicm = 0.0
      for j in range(monomers):
        p = i * monomers + j
        xi.append(x[p])
        yi.append(y[p])
        zi.append(z[p])
        xicm = xicm + x[p]
        yicm = yicm + y[p]
        zicm = zicm + z[p]
      xicm = xicm / monomers
      yicm = yicm / monomers
      zicm = zicm / monomers

      # store initial center-of-mass position
      com_zero.append((xicm, yicm, zicm))

  # close the file
  f.close()

  # initialize values
  g1 = 0.0
  g2 = 0.0
  g3 = 0.0

  # loop over chains
  for i in range(chains):

    # compute center-of-mass
    xcm = 0.0; ycm = 0.0; zcm = 0.0
    for j in range(monomers):
      p = i * monomers + j
      xcm = xcm + x[p]
      ycm = ycm + y[p]
      zcm = zcm + z[p]
      g1 = g1 + (x[p] - xi[p])**2 + (y[p] - yi[p])**2 + (z[p] - zi[p])**2
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers

    # compute g2
    for j in range(monomers):
      p = i * monomers + j
      g2 = g2 + (x[p] - xi[p] - xcm + com_zero[i][0])**2 + \
                (y[p] - yi[p] - ycm + com_zero[i][1])**2 + \
                (z[p] - zi[p] - zcm + com_zero[i][2])**2

    # compute g3
    g3 = g3 + (xcm - com_zero[i][0])**2 + (ycm - com_zero[i][1])**2 + (zcm - com_zero[i][2])**2

  # store data for the given file in the list
  time_g123.append((timeval, g1 / total_particles, g2 / total_particles, g3 / chains))

# write out the data
outfile = 'g123M' + str(monomers) + '.dat'
fout = open(outfile, 'w')
fout.write('#time (tau) g1 g2 g3 (sigma**2)\n')
for i in range(len(blockfiles)):
  fout.write('%d %.3f %.3f %.3f\n' % (time_g123[i][0], time_g123[i][1], time_g123[i][2], time_g123[i][3]))
fout.close()
sys.stdout.write(outfile + ' has been written to disk (time versus MSD).\n')

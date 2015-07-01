#!/usr/bin/env python

"""This Python script computes the average number of neighbors
   within Re of each chain as determined by the distance between
   the center of mass of each chain."""

# use sys for printing compatibility between Python 2 & 3
import sys
import numpy as np

# all systems have two hundred chains
chains = 200

# monomer number density (sigma**-3)
rho = 0.85

# intermediate output ('yes' or 'no')
verbose = 'yes'

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('*00000')

  # sort blockfiles list
  blockfiles.sort(lambda u, v: cmp(int(u[u.index('.') + 1:]), int(v[v.index('.') + 1:])))

  # work with subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 140000000
  iend = 440000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RGEQ_N200M400.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M400.0548100000', 'RGEQ_N200M400.0548190000']
 
else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

# initialize list
t_aveNeighbors = []

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
  L = (total_particles / rho)**(1.0 / 3.0)
  L_half = L / 2.0
  monomers = total_particles / chains

  # dictionary of N versus Rgsq and Resq
  N_Rgsq = {100:17.09, 200:29.73, 400:52.39, 800:85.62, 1600:146.23}
  N_Resq = {100:50.54, 200:84.66, 400:150.36, 800:245.69, 1600:388.67}

  # extract Re from dictionary
  Resq = N_Resq[monomers]

  # safety checks
  if(file == blockfiles[0]):
    monomers1st = monomers
    sys.stdout.write('chains = ' + str(chains) + ', monomers = ' + str(monomers) + '\n\n')
  if(monomers != monomers1st):
    sys.stderr.write('ERROR: monomers is not equal to monomers1st\n')
    sys.exit(1)
  if(verbose == 'yes' and file == blockfiles[0]):
      sys.stdout.write('Parameters found:\n')
      sys.stdout.write('  chains = ' + str(chains) + '\n' \
                       '  monomers = ' + str(monomers) + '\n' \
                       '  total_particles = ' + str(total_particles) + '\n' \
                       '  L = %.3f sigma\n  rho = %4.2f sigma**-3\n' % (L, rho))

  sys.stdout.write('Working on '+ file + '\n')

  # find line where time is found
  line = ''
  while not '{time ' in line:
    line = f.readline()
  sline = line.replace('}','')
  timeval = int(eval(sline.split()[1]))

  # check that time matches file extension (tau is 100 dt)
  exten = int(file[file.find('.0') + 1:]) / 100
  if(exten != timeval):
    sys.stdout.write('WARNING: Mismatch between time and file name ...\n')
    sys.stdout.write('WARNING: Using extension value of ' + str(exten) + \
                     ' instead of ' + str(timeval) + '\n\n') 
    timeval = exten

  # find line where coordinates begin
  line = ''
  while not '{particles ' in line:
    line = f.readline()

  # store the coordinates
  x = []; y = []; z = []
  for i in range(total_particles):
    r = f.readline().split()
    x.append(float(r[1]))
    y.append(float(r[2]))
    z.append(float(r[3]))

  # close the file
  f.close()

  com = []
  for i in range(chains):
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

    # shift chains such that the central box is between [0, L]
    ix = int(xcm / L)
    iy = int(ycm / L)
    iz = int(zcm / L)

    # correct in case of negative positions
    if(xcm < 0.0): ix = ix - 1
    if(ycm < 0.0): iy = iy - 1
    if(zcm < 0.0): iz = iz - 1

    # store shifted center-of-mass
    xcm = xcm - ix * L
    ycm = ycm - iy * L
    zcm = zcm - iz * L

    # store initial center-of-mass position
    com.append((xcm, ycm, zcm))

  count = 0
  # count neighbors
  for i in range(chains)[0:-1]:
    xi = com[i][0]
    yi = com[i][1]
    zi = com[i][2]
    for j in range(chains)[i+1:]:
      xij = xi - com[j][0]
      yij = yi - com[j][1]
      zij = zi - com[j][2]

      # apply MIC
      if(xij >  L_half): xij = xij - L;
      if(xij < -L_half): xij = xij + L;
      if(yij >  L_half): yij = yij - L;
      if(yij < -L_half): yij = yij + L;
      if(zij >  L_half): zij = zij - L;
      if(zij < -L_half): zij = zij + L;

      rijsq = xij**2 + yij**2 + zij**2
      if(rijsq < Resq): count = count + 1

  t_aveNeighbors.append((timeval, 2.0 * count / chains))

# compute average over all files
ct = map(lambda u: u[1], t_aveNeighbors)
ave_neighbors = np.mean(ct)
std_neighbors = np.std(ct)

# write data
outfile = 'ave_neighbors_' + str(monomers) + 'a.dat'
f = open(outfile, 'w')
f.write('# average number of neighbors within Re at t\n')
f.write('# <neighbors> = %10.3f, neighbors_std = %10.3f\n' % (ave_neighbors, std_neighbors))
f.write('# t (tau)  <neighbors>\n')
for i in range(len(t_aveNeighbors)):
  f.write('%10d %10.3f\n' % (t_aveNeighbors[i][0], t_aveNeighbors[i][1]))
f.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

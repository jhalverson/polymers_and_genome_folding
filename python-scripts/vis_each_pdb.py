"""This script loops over the specified block
   files and produces a PDB file for each chain
   up to show_chains sorted by zcm. Alternatively,
   one may also specify a list of chains. Wrapped
   coordinates are used. A PyMOL script is generated.
   The chains may be ouputted in either the monomer
   or center-of-mass particle representation.
"""

# use sys for printing compatibility between Python 2 & 3
import sys
import os
import random

# software (LAMMPS or else)
code = 'LAMMPS'

# all systems have two hundred chains
chains = 200

# density
rho = 0.85

# display intermediate output ('yes' or 'no')
verbose = 'yes'

# load output into PyMOL ('yes' or 'no')
pymol_load = 'no'

# PDB format string
pdb = '%6s%5d%5s%4s%6d%12.3f%8.3f%8.3f\n'

# scale factor for coordinates
sf = 2.0

# write box file ('yes' or 'no')
write_box = 'yes'

# how many of the chains should be shown
show_chains = chains

# track these chains only (make this empty to use zcm method)
track = []

# representation ('monomers' or 'com_particle')
rep = 'monomers'

# generation of file names: (1) glob or (2) loop or (3) manual
blockfiles_method = 3

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('ppa800a.*')

  # sort blockfiles list
  blockfiles.sort(lambda u, v: cmp(int(u[u.index('.', 5) + 1:]), int(v[v.index('.', 5) + 1:])))

  # work with subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 109000
  iend = 110000
  incr = 100

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('ppa800a.' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of block files: ' + str(len(blockfiles)) + '\n')

# seed the RNG
random.seed(10)

# create list of colors
colors = []
for i in range(chains):
  r = random.random()
  g = random.random()
  b = random.random()
  colors.append('[%4.3f, %4.3f, %4.3f]' % (r, g, b))

# set view
view = """
set_view (\
     1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000000,    0.000000000, -478.464263916,\
    60.532840729,   56.964488983,   57.608570099,\
   377.224884033,  579.703674316,    0.000000000 )
"""

if(code != 'LAMMPS'):
  # get value for monomers from file name
  prefix = blockfiles[0][0:blockfiles[0].index('.', blockfiles[0].index('R'))]
  monomers_filename = prefix[prefix.index('M') + 1:]
  if(monomers_filename[-1] != '0'):
    monomers_filename = monomers_filename[:-1]

# create empty PyMOL input file
fname = 'M.pml'
fpymol = open(fname, 'w')
fpymol.close()

# loop over block files
for file in blockfiles:

  f = open(file)

  if(code == 'LAMMPS'):

    sys.stdout.write('Working on '+ file + '\n')

    # find and store time
    line = ''
    while not 'ITEM: TIMESTEP' in line:
      line = f.readline()
    timeval = int(f.readline())

    # find and store total number of particles
    line = ''
    while not 'ITEM: NUMBER OF ATOMS' in line:
      line = f.readline()
    total_particles = int(f.readline())
    monomers = total_particles / chains

    # find and store box bounds
    line = ''
    while not 'ITEM: BOX BOUNDS' in line:
      line = f.readline()
    xmin, xmax = map(float, f.readline().split())
    # L = xmax - xmin
    L = (total_particles / rho)**(1.0 / 3.0)

    # find and store coordinates
    line = ''
    while not 'ITEM: ATOMS' in line:
      line = f.readline()
    coords = []
    for i in range(total_particles):
      sline = f.readline().split()
      coords.append((int(sline[0]), float(sline[1]), float(sline[2]), float(sline[3])))

    # sort the particles
    coords.sort(lambda u, v: cmp(u[0], v[0]))

    # store in linear lists
    x = []; y = []; z = []
    for i in range(total_particles):
      x.append(coords[i][1])
      y.append(coords[i][2])
      z.append(coords[i][3])

    # close the file
    f.close()

  else:

    # find line where total particles is found
    line = ''
    while not '{n_part ' in line:
      line = f.readline()

    sline = line.replace('}','')
    total_particles = int(sline.split()[1])
    L = (total_particles / rho)**(1.0 / 3.0)
    monomers = total_particles / chains

    if(file == blockfiles[0]):
      monomers1st = monomers
      if(monomers != int(monomers_filename)):
	sys.stdout.write('ERROR: monomers is not equal to monomers_filename\n\n')
      if(verbose == 'yes'):
	sys.stdout.write('Parameters found:\n')
	sys.stdout.write('  chains = ' + str(chains) + '\n' \
			 '  monomers = ' + str(monomers) + '\n' \
			 '  total_particles = ' + str(total_particles) + '\n' \
			 '  L = %g sigma\n  rho = %g sigma**-3\n\n' % (L, rho))
      if(len(track) != 0):
	sys.stdout.write('Tracking chains: ' + str(track[0]) + ', ' \
					     + str(track[1]) + ', ...\n\n')

    if(monomers != monomers1st or monomers != int(monomers_filename)):
      sys.stdout.write('ERROR: monomers is not equal to monomers1st\n')
      sys.exit(1)

    sys.stdout.write('Working on ' + file + ' ...\n')

    # find line where time is found
    line = ''
    while not '{time ' in line:
      line = f.readline()

    sline = line.replace('}','')
    timeval = int(eval(sline.split()[1]))

    # check that time matches file extension (tau is 100 dt)
    exten = int(file[file.find('.', file.index('R')) + 1:]) / 100
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
      sline = f.readline().split()
      x.append(float(sline[1]))
      y.append(float(sline[2]))
      z.append(float(sline[3]))
    f.close()

  # initialize list to store chain index and center-of-mass position
  com = []

  # loop over chains to wrap coordinates
  for i in range(chains):

    for j in range(monomers):
       p = i * monomers + j
       x[p] -= xmin
       y[p] -= xmin
       z[p] -= xmin

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

    # shift chain such that the central box is between [0, L]
    ix = int(xcm / L)
    iy = int(ycm / L)
    iz = int(zcm / L)
 
    # correct in case of negative positions
    if(xcm < 0.0): ix = ix - 1
    if(ycm < 0.0): iy = iy - 1
    if(zcm < 0.0): iz = iz - 1

    # wrap center of mass position
    xcm = xcm - ix * L
    ycm = ycm - iy * L
    zcm = zcm - iz * L

    com.append((i, xcm, ycm, zcm))
  
    # perform wrap
    for j in range(monomers):
      p = i * monomers + j
      x[p] = x[p] - ix * L
      y[p] = y[p] - iy * L
      z[p] = z[p] - iz * L

  # write out wrapped coords
  if 1:
    fout = open('coords.in', 'w')
    xminP = 1e6; yminP = 1e6; zminP = 1e6
    xmaxP = -1e6; ymaxP = -1e6; zmaxP = -1e6
    for i in range(len(x)):
      if(x[i] > xmaxP): xmaxP = x[i]
      if(y[i] > ymaxP): ymaxP = y[i]
      if(z[i] > zmaxP): zmaxP = z[i]
      if(x[i] < xminP): xminP = x[i]
      if(y[i] < yminP): yminP = y[i]
      if(z[i] < zminP): zminP = z[i]
    shft = 0.1
    print xminP, xmaxP, xmaxP - xminP + shft
    print yminP, ymaxP, ymaxP - yminP + shft
    print zminP, zmaxP, zmaxP - zminP + shft
    for i in range(len(x)):
      fout.write('%d %10.3f %10.3f %10.3f\n' % (i + 1, x[i] - xminP + shft, y[i] - yminP + shft, z[i] - zminP + shft))
    fout.close()

  # arrange com by zcm or track
  if(len(track) != 0):
    show_chains = len(track)
    track.sort()
    for i in range(len(track)):
      com[i] = com[track[i]]
    com = com[0:show_chains]
  else:
    # reverse sort the list and create sublist
    # com.sort(lambda u, v: cmp(u[3], v[3]), reverse = True)
    com = com[0:show_chains]

  # loop over chains and write a file for each
  for k in range(show_chains):
    # fname = 'M' + str(monomers) + 'T' + str(timeval) + 'C' + str(k) + '.pdb'
    fname = 'M' + str(monomers) + 'C' + str(k) + '.pdb'
    outfile = open(fname, 'w')
    i = com[k][0]
    if(rep == 'monomers'):
      # write monomers for each chain
      for j in range(monomers):
        p = i * monomers + j
        outfile.write(pdb % ('HETATM', j, 'C', 'ATM', i, sf * x[p], sf * y[p], sf * z[p]))
    else:
      # write center-of-mass particle for each chain 
      xcm = com[k][1]
      ycm = com[k][2]
      zcm = com[k][3]
      outfile.write(pdb % ('HETATM', 1, 'C', 'ATM', i, sf * xcm, sf * ycm, sf * zcm))
    outfile.close()
    if(verbose == 'yes'): sys.stdout.write('Wrote file ' + fname + '\n')

  # write PyMOL input file
  fname = 'M.pml'
  fpymol = open(fname, 'a')
  fpymol.write('reinitialize\n')
  fpymol.write('viewport 640,480\n')
  for k in range(show_chains):
    # flnm = 'M' + str(monomers) + 'T' + str(timeval) + 'C' + str(k) + '.pdb'
    flnm = 'M' + str(monomers) + 'C' + str(k) + '.pdb'
    fpymol.write('load ' + flnm + '\n')
    fpymol.write('set_color color' + str(k) + ', ' + colors[com[k][0]] + '\n')
    fpymol.write('color color' + str(k) + ', ' + flnm[0:len(flnm) - 4] + '\n')
    fpymol.write('show spheres, ' + flnm[0:len(flnm) - 4] + '\n')
  if(rep == 'com_particle'):
    fpymol.write('alter elem C, vdw=%.1f\n' % (monomers / 50))
    fpymol.write('rebuild\n')
  if(write_box == 'yes'):
    fpymol.write('load box.pdb\n')
    fpymol.write('color white, box\n')
    fpymol.write('show sticks, box\n')
  fpymol.write('reset\n')
# fpymol.write(view + '\n')
# fpymol.write('ray\n')
# fpymol.write('png M' + str(monomers) + 'T' + str(timeval) + '.png\n')
  fpymol.close()

# add quit statement
#fname = 'M' + str(monomers) + '.pml'
#fpymol = open(fname, 'a')
#fpymol.write('quit\n')
#fpymol.close()

sys.stdout.write('Wrote out PyMOL file ' + fname + '\n')

# write out box file
if(write_box == 'yes'):
  scaled_L = sf * L
  fbox = open('box.pdb', 'w')
  fbox.write(pdb % ('HETATM', 1, 'C', 'ATM', 1, 0.0, 0.0, 0.0))
  fbox.write(pdb % ('HETATM', 2, 'C', 'ATM', 2, 0.0, 0.0, scaled_L))
  fbox.write(pdb % ('HETATM', 3, 'C', 'ATM', 3, scaled_L, scaled_L, 0.0))
  fbox.write(pdb % ('HETATM', 4, 'C', 'ATM', 4, scaled_L, scaled_L, scaled_L))
  fbox.write(pdb % ('HETATM', 5, 'C', 'ATM', 5, 0.0, scaled_L, 0.0))
  fbox.write(pdb % ('HETATM', 6, 'C', 'ATM', 6, 0.0, scaled_L, scaled_L)) 
  fbox.write(pdb % ('HETATM', 7, 'C', 'ATM', 7, scaled_L, 0.0, 0.0))
  fbox.write(pdb % ('HETATM', 8, 'C', 'ATM', 8, scaled_L, 0.0, scaled_L))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 1, 2))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 3, 4))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 1, 5))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 1, 7))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 2, 6))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 2, 8))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 3, 5))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 4, 6))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 4, 8))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 5, 6))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 7, 8))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 3, 7))
  fbox.close()
  sys.stdout.write('box.pdb written to disk.' + '\n\n')

# view chains in PyMOL
if(len(blockfiles) == 1 and pymol_load == 'yes'):
  cmd = 'env pymol ' + fname
  os.system(cmd)

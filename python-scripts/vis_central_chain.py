"""
This script creates the PDB files necessary to animate the motion
of a central chain (with prescribed id) and its closest neighbors.
"""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os
import random

# id of central chain
id = 37

# how many neighbors should be shown (zero is valid)
neighbors = 3

# all systems have two hundred chains
chains = 200

# density
rho = 0.85

# PDB format string
pdb = '%6s%5d%5s%4s%6d%12.3f%8.3f%8.3f\n'

# scale factor for coordinates in PDB files
sf = 2.0

# translate chains by center-of-mass vector of central chain ('yes' or 'no')
translate = 'yes'

# write box file ('yes' or 'no')
write_box = 'no'

# representation ('monomers' or 'com_particle')
rep = 'monomers'

# display intermediate output ('yes' or 'no')
verbose = 'yes'

# load output into PyMOL ('yes' or 'no')
pymol_load = 'no'

# generation of file names: (2) loop or (3) manual
blockfiles_method = 3

if(blockfiles_method == 2):

  istr = 360800000
  iend = 362000000
  incr = 100000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RG_N200M1600.0361200000', 'RG_N200M1600.0361300000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n')
sys.stdout.write('############################################################\n')
sys.stdout.write('#                   Single chain animation                 #\n')
sys.stdout.write('############################################################\n')
sys.stdout.write('\n')

sys.stdout.write('Total number of block files: ' + str(len(blockfiles)) + '\n')

# seed the RNG
random.seed(10)

# create a pseudo-random list of colors
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
     0.000000000,    0.000000000, -292.506866455,\
    -1.439052582,   -3.374851227,    1.803752899,\
   230.614639282,  354.399078369,    0.000000000 )
"""

# get value for monomers from file name
prefix = blockfiles[0][0:blockfiles[0].index('.', blockfiles[0].index('R'))]
monomers_filename = prefix[prefix.index('M') + 1:]
if(monomers_filename[-1] != '0'):
  monomers_filename = monomers_filename[:-1]

# create empty PyMOL input file
fname = 'M' + monomers_filename + '.pml'
fpymol = open(fname, 'w')
fpymol.close()

# loop over block files
for file in blockfiles:

  f = open(file)

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
      sys.exit(1)
    if(verbose == 'yes'):
      sys.stdout.write('Parameters found:\n')
      sys.stdout.write('  chains = ' + str(chains) + '\n' \
                       '  monomers = ' + str(monomers) + '\n' \
                       '  total_particles = ' + str(total_particles) + '\n' \
                       '  L = %.3f sigma\n  rho = %4.2f sigma**-3\n' % (L, rho))

  if(monomers != monomers1st or monomers != int(monomers_filename)):
    sys.stdout.write('ERROR: monomers is not equal to monomers1st\n')
    sys.stdout.write('ERROR: or monomers is not equal to monomers_filename.\n')
    sys.exit(1)

  sys.stdout.write('\nWorking on ' + file + ' ...\n')

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

  # store the coordinates
  x = []; y = []; z = []
  for i in range(total_particles):
    sline = f.readline().split()
    x.append(float(sline[1]))
    y.append(float(sline[2]))
    z.append(float(sline[3]))
  f.close()

  # loop over chains to correct positions
  for i in range(chains):

    # compute center-of-mass
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

    # perform shift
    for j in range(monomers):
      p = i * monomers + j
      x[p] = x[p] - ix * L
      y[p] = y[p] - iy * L
      z[p] = z[p] - iz * L

    if(i == id):
      # store shifted center-of-mass
      xcm_central = xcm - ix * L
      ycm_central = ycm - iy * L
      zcm_central = zcm - iz * L
      # compute and store radius of gyration
      radsq = 0.0
      for j in range(monomers):
        p = i * monomers + j
        radsq = radsq + (x[p] - xcm_central)**2 \
                      + (y[p] - ycm_central)**2 \
                      + (z[p] - zcm_central)**2
      RGsq = radsq / monomers
      RG = math.sqrt(RGsq)

  if(neighbors != 0):
    # create 26 translations
    sys.stdout.write('\nCreating the 26 translations ... ')
    ict = 0
    p = total_particles
    for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          if(i**2 + j**2 + k**2 != 0):
            ict = ict + 1
            for m in range(total_particles):
              icentral = p - ict * total_particles
              if(0 > icentral > total_particles): sys.stdout.write('ERROR: icentral\n')
              x.append(x[icentral] + i * L)
              y.append(y[icentral] + j * L)
              z.append(z[icentral] + k * L)
              colors.append(colors[icentral])
              p = p + 1
    sys.stdout.write('done.\n')

    # initialize list of separations
    separations = []

    # loop over all chains excluding id
    sys.stdout.write('Computing separation distances ... ')
    for m in range(27 * chains):
      if(m != id):
        # compute center-of-mass
        xcm = 0.0
        ycm = 0.0
        zcm = 0.0
        for j in range(monomers):
          p = m * monomers + j
          xcm = xcm + x[p]
          ycm = ycm + y[p]
          zcm = zcm + z[p]
        xcm = xcm / monomers
        ycm = ycm / monomers
        zcm = zcm / monomers

        rsq = (xcm_central - xcm)**2 + (ycm_central - ycm)**2 + (zcm_central - zcm)**2
        separations.append((m, rsq))
    sys.stdout.write('done.\n\n')

    # sort the list of separations by rsq
    separations.sort(lambda u, v: cmp(u[1], v[1]))

    # build string of nearest neighbors
    nn = ''
    for n in range(neighbors):
      rij = math.sqrt(separations[n][1])
      rijstr = '%.1f' % rij
      nn = nn + str(separations[n][0]) + '(' + rijstr + ' sigma) '
    sys.stdout.write('Nearest neighbors of ' + str(id) + ': ' + nn + '\n')

  # write PDB file of central chain
  fname = 'M' + str(monomers) + 'T' + str(timeval) + 'ID' + str(id) + '.pdb'
  outfile = open(fname, 'w')
  for j in range(monomers):
    p = id * monomers + j
    if(translate == 'yes'):
      xtns = x[p] - xcm_central
      ytns = y[p] - ycm_central
      ztns = z[p] - zcm_central
      outfile.write(pdb % ('HETATM', j, 'C', 'ATM', id, sf * xtns, sf * ytns, sf * ztns))
    else:
      outfile.write(pdb % ('HETATM', j, 'C', 'ATM', id, sf * x[p], sf * y[p], sf * z[p]))
  outfile.close()
  sys.stdout.write(fname + ' for id ' + str(id) + ' has been written to disk.\n')

  if(neighbors != 0):
    sys.stdout.write('  Nearest neighbors may have an index between 0 and ' + str(27 * chains - 1) + '\n')

  # write out neighbors
  for k in range(neighbors):
    # neighbor id
    nid = separations[k][0]
    # write PDB file
    fname = 'M' + str(monomers) + 'T' + str(timeval) + 'ID' + str(id) + 'NN' + str(k) + '.pdb'
    outfile = open(fname, 'w')
    for j in range(monomers):
      p = nid * monomers + j
      if(translate == 'yes'):
        xtns = x[p] - xcm_central
        ytns = y[p] - ycm_central
        ztns = z[p] - zcm_central
        outfile.write(pdb % ('HETATM', j, 'C', 'ATM', nid, sf * xtns, sf * ytns, sf * ztns))
      else:
        outfile.write(pdb % ('HETATM', j, 'C', 'ATM', nid, sf * x[p], sf * y[p], sf * z[p]))
    outfile.close()
    sys.stdout.write('  ' + fname + ' for nearest neighbor ' + str(nid) + ' of ' \
                     + str(id) + ' has been written to disk.\n')

  # write PyMOL input file
  fname = 'M' + str(monomers) + '.pml'
  fpymol = open(fname, 'a')
  fpymol.write('reinitialize\n')
  fpymol.write('viewport 640,480\n')
  flnm = 'M' + str(monomers) + 'T' + str(timeval) + 'ID' + str(id) + '.pdb'
  fpymol.write('load ' + flnm + '\n')
  fpymol.write('set_color color0, ' + colors[id] + '\n')
  fpymol.write('color color0, ' + flnm[0:len(flnm) - 4] + '\n')
  fpymol.write('show spheres, ' + flnm[0:len(flnm) - 4] + '\n')
  for k in range(neighbors):
    nid = separations[k][0]
    flnm = 'M' + str(monomers) + 'T' + str(timeval) + 'ID' + str(id) + 'NN' + str(k) + '.pdb'
    fpymol.write('load ' + flnm + '\n')
    fpymol.write('set_color color' + str(k + 1) + ', ' + colors[nid] + '\n')
    fpymol.write('color color' + str(k + 1) + ', ' + flnm[0:len(flnm) - 4] + '\n')
    fpymol.write('show spheres, ' + flnm[0:len(flnm) - 4] + '\n')
  if(rep == 'com_particle'):
    fpymol.write('alter elem C, vdw=%.1f\n' % (monomers / 100))
    fpymol.write('rebuild\n')
  if(write_box == 'yes'):
    fpymol.write('load box.pdb\n')
    fpymol.write('color white, box\n')
    fpymol.write('show sticks, box\n')
  fpymol.write(view + '\n')
  fpymol.write('ray\n')
  fpymol.write('png M' + str(monomers) + 'T' + str(timeval) + 'RG%.1f.png\n' % (RG))
  fpymol.close()

# add quit statement
fname = 'M' + str(monomers) + '.pml'
fpymol = open(fname, 'a')
fpymol.write('quit\n')
fpymol.close()

sys.stdout.write('\n' + fname + ' written to disk.\n')

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
  sys.stdout.write('box.pdb written to disk.\n\n')

# view chains in PyMOL
if(len(blockfiles) == 1 and pymol_load == 'yes'):
  cmd = 'env pymol ' + fname
  os.system(cmd)

"""This script computes the radius of gyration of each chain
   and writes the n largest chains and their nearest neighbors
   to a PDB file for visualization."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os

# how many of the largest chains should be considered
nlargest = 1

# how many of the nearest neighbors of the largest chains should be shown
nnshow = 2

# all systems have two hundred chains
chains = 200

# density
rho = 0.85

# display intermediate output ('yes' or 'no')
verbose = 'yes'

# load output into PyMol ('yes' or 'no')
pymol_load = 'yes'

# generation of file names: (1) glob, (2) loop, (3) manual
blockfiles_method = 3

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('*000000')

  # sort blockfiles list
  blockfiles.sort(lambda u, v: cmp(int(u[u.index('.') + 1:]), int(v[v.index('.') + 1:])))

elif(blockfiles_method == 2):

  istr = 360800000
  iend = 362000000
  incr = 100000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RG_N200M1600.0361200000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of block files: ' + str(len(blockfiles)) + '\n')

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
    if(verbose == 'yes'):
      sys.stdout.write('Parameters found:\n')
      sys.stdout.write('  chains = ' + str(chains) + '\n' \
                       '  monomers = ' + str(monomers) + '\n' \
                       '  total_particles = ' + str(total_particles) + '\n' \
                       '  L = %g sigma\n  rho = %g sigma**-3\n\n' % (L, rho))

  if(monomers != monomers1st):
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

  # initialize list of pairs (i, RG)
  pair = []

  # loop over chains to find rcm and shift coordinates
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

    RGsq = radsq / monomers
    RG = math.sqrt(RGsq)
    pair.append((i, RG))

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

  # reverse sort the pairs by RG
  pair.sort(lambda u, v: cmp(u[1], v[1]), reverse = True)

  # write out i and RG of the n largest chains
  sys.stdout.write('The ' + str(nlargest) + ' largest chains: \n') 
  for i in range(nlargest):
    sys.stdout.write('  index = %d, RG = %g sigma \n' % (pair[i][0], pair[i][1]))

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
            p = p + 1
  sys.stdout.write('done.\n\n')

  # loop over n largest chains and find nearest neighbors
  for i in range(nlargest):

    # store index of large chain
    ilarge = pair[i][0]

    # compute center-of-mass
    xcm = 0.0; ycm = 0.0; zcm = 0.0
    for j in range(monomers):
      p = ilarge * monomers + j
      xcm = xcm + x[p]
      ycm = ycm + y[p]
      zcm = zcm + z[p]
    xcm_large = xcm / monomers
    ycm_large = ycm / monomers
    zcm_large = zcm / monomers

    # initialize nearest neighbor list
    nnpair = []

    # loop over all chains excluding ilarge
    for m in range(27 * chains):
      if(m != ilarge):

        # compute center-of-mass
        xcm = 0.0; ycm = 0.0; zcm = 0.0
        for j in range(monomers):
          p = m * monomers + j
          xcm = xcm + x[p]
          ycm = ycm + y[p]
          zcm = zcm + z[p]
        xcm = xcm / monomers
        ycm = ycm / monomers
        zcm = zcm / monomers

        rsq = (xcm_large - xcm)**2 + (ycm_large - ycm)**2 + (zcm_large - zcm)**2
        nnpair.append((m, rsq))

    # sort the nearest neighbor pairs by RG
    sys.stdout.write('Sorting nearest neighbors by distance ... ')
    nnpair.sort(lambda u, v: cmp(u[1], v[1]))
    sys.stdout.write('done.\n\n')

    # build string of nearest neighbors
    nn = ''
    for n in range(nnshow):
      rij = math.sqrt(nnpair[n][1])
      rijstr = '%g' % rij
      nn = nn + str(nnpair[n][0]) + '(' + rijstr + ' sigma) '

    sys.stdout.write('Nearest neighbors of large chain ' \
                     + str(ilarge) + ': ' + nn + '\n')

    # PDB format string
    pdb = '%6s%5d%5s%4s%6d%12.3f%8.3f%8.3f\n'

    # scale factor for coordinates
    sf = 2.0

    # initialize list of PDB file names
    flnm_pdb = []
 
    # write PDB file
    fname = 'M' + str(monomers) + 'T' + str(timeval) + 'L' + str(i) + '.pdb'
    flnm_pdb.append(fname)
    outfile = open(fname, 'w')

    # write out central chain
    for j in range(monomers):
      p = ilarge * monomers + j
      outfile.write(pdb % ('HETATM', ilarge, 'C', 'ATM', j, sf * x[p], sf * y[p], sf * z[p]))
    outfile.close()
    sys.stdout.write(fname + ' for large chain ' + str(ilarge) + ' has been written to disk.\n')
    sys.stdout.write('  Nearest neighbors may have an index between 0 and ' + str(27 * chains - 1) + '\n')

    # write out nearest neighbors
    for k in range(nnshow):

      # write PDB file
      fname = 'M' + str(monomers) + 'T' + str(timeval) + 'L' + str(i) + 'NN' + str(k) + '.pdb'
      flnm_pdb.append(fname)
      outfile = open(fname, 'w')

      for j in range(monomers):
        p = nnpair[k][0] * monomers + j
        outfile.write(pdb % ('HETATM', nnpair[k][0], 'C', 'ATM', j, sf * x[p], sf * y[p], sf * z[p]))
    
      outfile.close()
      sys.stdout.write('  ' + fname + ' for nearest neighbor ' + str(nnpair[k][0]) + ' of large chain ' \
                       + str(ilarge) + ' has been written to disk.\n')

    # list of colors for Pymol
    colors = ['red', 'green', 'yellow', 'cyan', 'orange', 'white', 'gray', 'magenta']

    if(nlargest == 1 and len(colors) >= 1 + nnshow):

      # write PyMol input file
      fname = 'M' + str(monomers) + 'T' + str(timeval) + 'L' + str(i) + '.pml'
      fpymol = open(fname, 'w')
      fpymol.write('reinitialize\n')
      fpymol.write('viewport 640,480\n')
      for k, flnm in enumerate(flnm_pdb):
        fpymol.write('load ' + flnm + '\n')
        fpymol.write('color ' + colors[k] + ', ' +  flnm[0:len(flnm) - 4] + '\n')
        fpymol.write('show spheres, ' + flnm[0:len(flnm) - 4] + '\n')
      fpymol.write('reset\n')
      fpymol.close()
      sys.stdout.write('\nWrote out PyMol file ' + fname + '\n')

      if(len(blockfiles) == 1 and pymol_load == 'yes'):
        # view chains in PyMol
        cmd = 'env pymol ' + fname
        os.system(cmd)

  sys.stdout.write('\nDone with file ' + file + '\n\n')

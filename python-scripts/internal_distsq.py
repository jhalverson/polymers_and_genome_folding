"""This Python script computes the probability of monomer i
   interacting with monomer j of the same chain."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os

# all systems have two hundred chains
chains = 200

# ESPResSo or LAMMPS
code = 'LAMMPS'

# rings or linear
topology = 'rings'

# intermediate output ('yes' or 'no')
verbose = 'yes'

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
# blockfiles = glob.glob('*0000')
  blockfiles = glob.glob('rings800a.*')

  # sort blockfiles list
  blockfiles.sort(lambda u, v: cmp(int(u[u.index('.', 5) + 1:]), int(v[v.index('.', 5) + 1:])))

  # work with subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M400.0450870000']
 
else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write(' ' + code + '\n')
sys.stdout.write('  Re**2 will be computed with topology = ' + topology + '\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

r2 = []
norm = []
for i in xrange(400+1):
  r2.append(0)
  norm.append(0)
sys.stdout.write('monomers_half = 400?\n\n')

# loop over files
for file in blockfiles:

  if(code == 'LAMMPS'):

    # open file
    f = open(file)

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
    monomers_half = monomers / 2

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

    # open file
    f = open(file)

    # find line where total particles is found
    line = ''
    while not '{n_part ' in line:
      line = f.readline()
    sline = line.replace('}','')
    total_particles = int(sline.split()[1])
    monomers = total_particles / chains
    monomers_half = monomers / 2

    # safety checks
    if(file == blockfiles[0]):
      monomers1st = monomers
      sys.stdout.write('chains = ' + str(chains) + ', monomers = ' + str(monomers) + '\n\n')
    if(monomers != monomers1st):
      sys.stderr.write('ERROR: monomers is not equal to monomers1st\n')
      sys.exit(1)

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
      sys.stderr.write('WARNING: Mismatch between time and file name ...\n')
      sys.stderr.write('WARNING: Using extension value of ' + str(exten) + \
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

    # close the file
    f.close()

# the two paths remerge

  def nbonds_separated(a_, b_):
    """Assumes particle indices are between 0 and monomers - 1"""
    if(a_ <= monomers_half and b_ <= monomers_half): return abs(a_ - b_)
    if(a_ >= monomers_half and b_ >= monomers_half): return abs(a_ - b_)
 
    if(a_ >= monomers_half and b_ <= monomers_half):
      sep1 = (monomers_half - b_) + (a_ - monomers_half)
      sep2 = b_ + (monomers - a_)
      return min(sep1, sep2)

    if(a_ <= monomers_half and b_ >= monomers_half):
      sep1 = (monomers_half - a_) + (b_ - monomers_half)
      sep2 = a_ + (monomers - b_)
      return min(sep1, sep2)

  """      
    monomers = 8
    monomers_half = 4
    print nbonds_separated(0, 1), nbonds_separated(0, 4), nbonds_separated(2, 4)
    print nbonds_separated(1, 0), nbonds_separated(4, 2), nbonds_separated(3, 0)
    print nbonds_separated(4, 5), nbonds_separated(7, 5), nbonds_separated(5, 7)
    print nbonds_separated(0, 7), nbonds_separated(3, 5), nbonds_separated(1, 6)
    print nbonds_separated(7, 0), nbonds_separated(5, 3), nbonds_separated(6, 1)
    print nbonds_separated(0, 4), nbonds_separated(4, 0)
    sys.exit()
  """

  for k in xrange(chains):
    for i in xrange(0, monomers - 1):
      p = k * monomers + i
      xi = x[p]
      yi = y[p]
      zi = z[p]
      for j in xrange(i + 1, monomers):
        p = k * monomers + j
        dsq = (xi - x[p])**2 + (yi - y[p])**2 + (zi - z[p])**2
        nij = nbonds_separated(i, j)
        r2[nij] = r2[nij] + dsq
        norm[nij] += 1

outfile = 'r2_' + str(monomers) + 'a.dat'
fout = open(outfile, 'w')
fout.write('# n  count  r2/n\n')
for i in xrange(1, monomers_half+1):
  fout.write('%d %d %8.6f\n' % (i, norm[i], r2[i]/norm[i]))
fout.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

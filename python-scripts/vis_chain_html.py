"""This script creates a PyMOL input file and an HTML
   file to display images of individual chains. The
   list of chains to display is appended by the chains
   with the smallest and largest radius of gyration.
"""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os
import random

# all systems have two hundred chains
chains = 200

# density
rho = 0.85

# display intermediate output ('yes' or 'no')
verbose = 'yes'

# load output into PyMol ('yes' or 'no')
pymol_load = 'yes'

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('100/*000')

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

  blockfiles = ['RGEQ_N200M800a.012400000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

# PDB format string
pdb = '%6s%5d%5s%4s%6d%12.3f%8.3f%8.3f\n'

# scale factor for coordinates
sf = 2.0

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
     0.000000000,   -0.000000045, -356.651336670,\
    -0.240227178,    3.604017973,    0.000076294,\
  -5249.776855469, 5963.079101562,    0.000000000 )
"""

# which chains should be written: (1) all, (2) every nth, (3) random or (4) only a few
chains2file = 3

if(chains2file == 1):
  # all chains will be written
  indices = range(chains)
elif(chains2file == 2):
  # every nth chain
  nth = 20
  indices = range(0, chains, nth)
elif(chains2file == 3):
  # n random chains
  n = 10
  indices = random.sample(xrange(chains), n)
elif(chains2file == 4):
  # only these chains
  indices = [0, 37, 40]
else:
  sys.stdout.write('Value of chains2file is not valid.\n\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of block files: ' + str(len(blockfiles)) + '\n')

# loop over block files
for a, file in enumerate(blockfiles):

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

  # initialize list (i, RG, xcm, ycm, zcm)
  props = []

  # loop over chains to find RG and center-of-mass
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

    # sum the squares
    radsq = 0.0
    for j in range(monomers):
      p = i * monomers + j
      radsq = radsq + (x[p] - xcm)**2 + (y[p] - ycm)**2 + (z[p] - zcm)**2

    # compute radius of gyration
    RGsq = radsq / monomers
    RG = math.sqrt(RGsq)
    props.append((i, RG, xcm, ycm, zcm))

  # sort props by RG
  props.sort(lambda u, v: cmp(u[1], v[1]))

  # write out properties of the smallest and largest chains
  sys.stdout.write('  Smallest: %d, RG = %.1f sigma\n' % (props[0][0], props[0][1])) 
  sys.stdout.write('  Largest:  %d, RG = %.1f sigma\n' % (props[-1][0], props[-1][1]))

  # add smallest to beginning and largest to end
  indices.insert(0, props[0][0])
  indices.append(props[-1][0])

  # resort by index for easy indexing
  props.sort(lambda u, v: cmp(u[0], v[0]))

  # initial list of PDB file names
  flnm_pdb = []

  # write a PDB file for each chain in indices
  for c in indices:
    RG  = props[c][1]
    xcm = props[c][2]
    ycm = props[c][3]
    zcm = props[c][4]
    fname = 'M' + str(monomers) + 'T' + str(timeval) + 'ID%dRG%.1f.pdb' % (c, RG)
    flnm_pdb.append(fname)
    outfile = open(fname, 'w')
    for j in range(monomers):
      p = c * monomers + j
      outfile.write(pdb % ('HETATM', j, 'C', 'ATM', c, \
                    sf * (x[p] - xcm), sf * (y[p] - ycm), sf * (z[p] - zcm)))
    outfile.close()

  # write PyMol input file
  fname = 'M' + str(monomers) + 'T' + str(timeval) + 'F' + str(a) + '.pml'
  fpymol = open(fname, 'w')
  fpymol.write('reinitialize\n')
  fpymol.write('viewport 480,480\n')
  for i, c in enumerate(indices):
    flnm = flnm_pdb[i]
    fpymol.write('set ray_trace_mode,2\n')
    fpymol.write('cmd.bg_color(\'white\')\n')
    fpymol.write('load ' + flnm + '\n')
    fpymol.write('show spheres, ' + flnm[0:len(flnm) - 4] + '\n')
    fpymol.write(view + '\n')
    fpymol.write('ray\n')
    fpymol.write('png ' + flnm[0:len(flnm) - 4] + 'z.png\n')
    fpymol.write('rotate y, 90\n')
    fpymol.write('ray\n')
    fpymol.write('png ' + flnm[0:len(flnm) - 4] + 'y.png\n')
    fpymol.write('rotate y, -90\n')
    fpymol.write('rotate x, 90\n')
    fpymol.write('ray\n')
    fpymol.write('png ' + flnm[0:len(flnm) - 4] + 'x.png\n')
    fpymol.write('reinitialize\n')
  fpymol.write('quit\n')
  fpymol.close()
  sys.stdout.write('  Wrote out PyMol file ' + fname + '\n')

  # view chains in PyMol
  if(len(blockfiles) == 1 and pymol_load == 'yes'):
    cmd = 'env pymol ' + fname
    os.system(cmd)

  # write HTML file
  fname = 'M' + str(monomers) + 'T' + str(timeval) + 'F' + str(a) + '.html'
  fhtml = open(fname, 'w')
  fhtml.write('<html>\n')
  fhtml.write('<head><title>' + str(monomers) + 'a</title></head>\n')
  fhtml.write('<body><center><table cellspacing="2" cellpadding="2">\n')
  for i, c in enumerate(indices):
    flnm = flnm_pdb[i]
    flnm = flnm[0:len(flnm) - 4]
    fhtml.write('<tr>\n')
    fhtml.write('<td><img src="' + flnm + 'z.png"></td>\n')
    fhtml.write('<td><img src="' + flnm + 'y.png"></td>\n')
    fhtml.write('<td><img src="' + flnm + 'x.png"></td>\n')
    fhtml.write('</tr>\n') 
  fhtml.write('</table></center></body></html>\n')
  fhtml.close()
  sys.stdout.write('  Wrote out HTML file ' + fname + '\n')

  sys.stdout.write('\nDone with file ' + file + '\n\n')

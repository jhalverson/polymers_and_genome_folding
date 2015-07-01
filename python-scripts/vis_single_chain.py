"""This Python script extracts a single chain from the system
   and creates a PNG image of the chain for each configuration
   file."""

import sys
import pypolymer

# which chain to make animation [0, chains - 1]
chain_id = 1

# time step (tau)
dt = 0.01

# ghost effect (0 equals no ghosts)
ghost_frames = 0

# COM fixed
com_fixed = True

# scale factor
sf = 2.0

# dictionary for Rg2
Rg2_ave = {100:17.16, 200:30.77, 400:52.92, 800:87.61, 1600:145.61}
#Rg2_ave = {100:43.35, 200:88.93, 400:180.84, 800:359.14}

# PDB format string
fmt = '%6s%5d%5s%4s%6d%12.3f%8.3f%8.3f\n'

# chain color
color0 = '[0.8, 0.8, 0.8]'

# set view
view = """
set_view (\
     1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000000,    0.000000000, -334.695098877,\
     0.245868564,   -2.458981276,    0.000007629,\
  -219.102279663,  888.492370605,    0.000000000 )
"""

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 2

if(blockfiles_method == 1):

  prefix = '../stripped/'
  existingfiles = glob.glob(prefix + '*.pos')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]

elif(blockfiles_method == 2):

  prefix = '../stripped/'
  istr = 472000000
  iend = 740000000
  incr = 1000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append(prefix + str(i) + '.pos')

elif(blockfiles_method == 3):

  blockfiles = ['../../']

else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit()

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                   Single chain                  \n')
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

data = []
for file in blockfiles:

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

  """
  total_particles = 320000
  M_rings = 0
  M_linear = 400
  monomers = 800
  monomers_half = 400
  """

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

  if(file == blockfiles[0]):
    xcm_sys_zero = 0.0
    ycm_sys_zero = 0.0
    zcm_sys_zero = 0.0
    for i in range(M_rings + M_linear):
      for j in range(monomers):
	p = i * monomers + j
	xcm_sys_zero = xcm_sys_zero + x[p]
	ycm_sys_zero = ycm_sys_zero + y[p]
	zcm_sys_zero = zcm_sys_zero + z[p]
    xcm_sys_zero = xcm_sys_zero / total_particles
    ycm_sys_zero = ycm_sys_zero / total_particles
    zcm_sys_zero = zcm_sys_zero / total_particles

  # compute COM of system
  xcm_sys = 0.0
  ycm_sys = 0.0
  zcm_sys = 0.0
  for i in range(M_rings + M_linear):
    for j in range(monomers):
      p = i * monomers + j
      xcm_sys = xcm_sys + x[p]
      ycm_sys = ycm_sys + y[p]
      zcm_sys = zcm_sys + z[p]
  xcm_sys = xcm_sys / total_particles
  ycm_sys = ycm_sys / total_particles
  zcm_sys = zcm_sys / total_particles

  # correct positions
  for i in range(M_rings + M_linear):
    for j in range(monomers):
      p = i * monomers + j
      x[p] = x[p] - (xcm_sys - xcm_sys_zero)
      y[p] = y[p] - (ycm_sys - ycm_sys_zero)
      z[p] = z[p] - (zcm_sys - zcm_sys_zero)

  xc = []
  yc = []
  zc = []
  for j in range(monomers):
    xc.append(x[chain_id * monomers + j])
    yc.append(y[chain_id * monomers + j])
    zc.append(z[chain_id * monomers + j])

  # compute center-of-mass
  xcm = 0.0
  ycm = 0.0
  zcm = 0.0
  for j in range(monomers):
    xcm = xcm + xc[j]
    ycm = ycm + yc[j]
    zcm = zcm + zc[j]
  xcm = xcm / monomers
  ycm = ycm / monomers
  zcm = zcm / monomers

  if(file == blockfiles[0]):
    xcm_zero = xcm
    ycm_zero = ycm
    zcm_zero = zcm

  # compute MSD
  msd = (xcm - xcm_zero)**2 + (ycm - ycm_zero)**2 + (zcm - zcm_zero)**2

  if(com_fixed):
    for j in range(monomers):
      xc[j] = xc[j] - xcm
      yc[j] = yc[j] - ycm
      zc[j] = zc[j] - zcm
    xcm = 0.0
    ycm = 0.0
    zcm = 0.0

  # compute Rg2
  Rg2 = 0.0
  for j in range(monomers):
    Rg2 = Rg2 + (xc[j] - xcm)**2 + (yc[j] - ycm)**2 + (zc[j] - zcm)**2
  Rg2 = Rg2 / monomers

  Rg_ratio = (Rg2 / Rg2_ave[monomers])**0.5
  g3_ratio = (msd / Rg2_ave[monomers])**0.5
  outfile = str(timestep) + '_' + str(chain_id) + '.pdb'
  data.append((timestep, Rg_ratio, g3_ratio, outfile))

  f = open(outfile, 'w')
  for j in range(monomers):
    f.write(fmt % ('HETATM', j + 1, 'X', 'PRT', 1, sf * xc[j], sf * yc[j], sf * zc[j]))
  f.close()

# write PyMOL input file
# opacity 1 is invisible while 0 is solid
fname = 'M.pml'
fpymol = open(fname, 'w')
for i, d in enumerate(data):
  fpymol.write('reinitialize\n')
  fpymol.write('viewport 640,480\n')
  fpymol.write('set_color color0, ' + color0 + '\n')
  j = min(i, ghost_frames)
  while(i - j >= 0 and j >= 0):
    flnm = data[i - j][3]
    fpymol.write('load ' + flnm + '\n')
    fpymol.write('color color0, ' + flnm[:-4] + '\n')
    fpymol.write('show spheres, ' + flnm[:-4] + '\n')
    fpymol.write('hide lines, ' + flnm[:-4] + '\n')
    if(ghost_frames != 0 and j != 0):
      trans = '%.3f' % (float(j) / (ghost_frames + 1))
      fpymol.write('set sphere_transparency=' + trans + ', ' + flnm[:-4] + '\n')
    j = j - 1
  fpymol.write('reset\n')
  fpymol.write(view + '\n')
  fpymol.write('ray\n')
  fpymol.write('png M' + str(monomers) + 'T' + str(d[0]) + '.png\n')
fpymol.write('quit\n')
fpymol.close()

# time written in tau
f = open('data.out', 'w')
for d in data:
  f.write('%d %.2f %.2f %s\n' % (d[0] * dt, d[1], d[2], d[3]))
f.close()

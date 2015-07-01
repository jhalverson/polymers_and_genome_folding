"""This script computes g1, g2 and g3; it tracks center-of-mass
   drift, allows one to average over multiple time origins and
   can output values uniformly spaced on a logarithmic scale."""

# use sys for printing compatibility between Python 2 & 3
import sys
import os
import time
import math
import glob
import pypolymer

# job case
case = ''

# inner monomers
inner_monomers = 5

# displacements for time origins (range may be used)
origins = []
origins = map(int, origins)

# generation of file names: (1) glob, (2) loop, (3) manual or (4) log scale
blockfiles_method = 4

if(blockfiles_method == 1):

  prefix = '../../rings200_26.'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]
  if 1:
    maxfiles = 25
    numfiles = len(blockfiles)
    if(numfiles > maxfiles): blockfiles = blockfiles[::int(numfiles / float(maxfiles))]

elif(blockfiles_method == 2):

  prefix = '../../rings200_26.'
  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append(prefix + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['file1', 'file2']

elif(blockfiles_method == 4):

  # prefix
  prefix = '../../rings200_26.'

  # maximum number of files (script may use less)
  n = 200

  # exponent of first and last decade (time steps)
  start_exp = 3
  end_exp = 10
 
  step_exp = float(end_exp - start_exp) / (n - 1)
  uniform_time = [10**(start_exp + step_exp * u) for u in range(n)]
  existingfiles = glob.glob(prefix + '*')
  if(len(existingfiles) == 0): sys.stderr.write('ERROR: No files found.\n'); sys.exit(1)
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles = [int(u[u.rindex('.') + 1:]) for u in existingfiles]
  existingfiles.sort()
  uniform_time = [pypolymer.find_closest_existing_file(existingfiles, u) for u in uniform_time]
  uniform_time = pypolymer.unique(uniform_time)
  blockfiles = [prefix + str(u) for u in uniform_time]

else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit()

if(len(blockfiles) == 0):
  sys.stderr.write('ERROR: Blockfiles is empty.\n')
  sys.exit()

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                        g123                     \n')
sys.stdout.write('=================================================\n')
sys.stdout.write('\n\nTotal number of files in first set: ' + str(len(blockfiles)) + '\n')
if(len(blockfiles) > 2):
  sys.stdout.write('  blockfiles[0] = ' + blockfiles[0] + '\n')
  sys.stdout.write('  blockfiles[1] = ' + blockfiles[1] + '\n')
  sys.stdout.write('  blockfiles[2] = ' + blockfiles[2] + '\n')
  sys.stdout.write('  ...' + '\n')
  sys.stdout.write('  blockfiles[-3] = ' + blockfiles[-3] + '\n')
  sys.stdout.write('  blockfiles[-2] = ' + blockfiles[-2] + '\n')
  sys.stdout.write('  blockfiles[-1] = ' + blockfiles[-1] + '\n\n')

# create list of file sets
fset = [blockfiles]

if(len(origins) != 0):
  # store existing files
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)

  # build file sets for the different time origins
  for i in range(len(origins)):
    fs = []
    disp = origins[i]
    for file in blockfiles:
      suffix = int(file[file.rindex('.') + 1:])
      fname = prefix + str(suffix + disp)
      fs.append(fname)
    fs = list(set(fs) & set(existingfiles))
    fs.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
    fset.append(fs)

sys.stdout.write('Number of file sets = ' + str(len(fset)) + '\n')
for i in range(len(fset)):
  sys.stdout.write('  number of files in set ' + str(i) + ' = ' + str(len(fset[i])) + '\n')

# initialize dictionaries for g123
g1ave_rings = {}
g2ave_rings = {}
g3ave_rings = {}
g1ave_linear = {}
g1ave_linear_inner = {}
g2ave_linear = {}
g3ave_linear = {}
cntr = {}

# for each file set
for files in fset:

  xi = []; yi = []; zi = []; com_zero = []
  sys.stdout.write('\n\nStarting new file set ...\n')

  # for each file in the file set
  for file in files:

    fmanager = pypolymer.read_blockfile(file)
    timestep = fmanager[0]
    total_particles = fmanager[1]
    M_rings = fmanager[2]; M_linear = fmanager[3]
    monomers = fmanager[4]; monomers_half = fmanager[5]
    x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

    if(files == fset[0] and file == files[0]):
      total_particles1st = total_particles
      sys.stdout.write('System parameters:\n')
      sys.stdout.write('  M_rings = ' + str(M_rings) + '\n')
      sys.stdout.write('  M_linear = ' + str(M_linear) + '\n')
      sys.stdout.write('  N = ' + str(monomers) + '\n')
      sys.stdout.write('  N/2 = ' + str(monomers_half) + '\n\n')
    if(total_particles1st != total_particles):
      sys.stderr.write('Change in total_particles (' + str(total_particles) + ').\n')
      sys.exit()

    sys.stdout.write('Working on '+ file + '\n')

    # use first file of each file set to get initial values
    if(file == files[0]):
      # store reference time
      time_zero = timestep
      if(timestep != 0):
	sys.stdout.write('WARNING: time is not 0 for initial file (' + str(timestep) +  ').\n\n')
      # loop over chains
      xcm_sys_zero = 0.0; ycm_sys_zero = 0.0; zcm_sys_zero = 0.0
      for i in range(M_rings + M_linear):
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
	  xcm_sys_zero = xcm_sys_zero + x[p]
	  ycm_sys_zero = ycm_sys_zero + y[p]
	  zcm_sys_zero = zcm_sys_zero + z[p]
	xicm = xicm / monomers
	yicm = yicm / monomers
	zicm = zicm / monomers

	# store initial center-of-mass position
	com_zero.append((xicm, yicm, zicm))

      xcm_sys_zero = xcm_sys_zero / total_particles
      ycm_sys_zero = ycm_sys_zero / total_particles
      zcm_sys_zero = zcm_sys_zero / total_particles
      sys.stdout.write('r_COM(t=' + str(timestep) + '): %8.2f %8.2f %8.2f\n' % \
		      (xcm_sys_zero, ycm_sys_zero, zcm_sys_zero))

    # compute COM of system
    xcm_sys = 0.0; ycm_sys = 0.0; zcm_sys = 0.0
    for i in range(M_rings + M_linear):
      for j in range(monomers):
	p = i * monomers + j
	xcm_sys = xcm_sys + x[p]
	ycm_sys = ycm_sys + y[p]
	zcm_sys = zcm_sys + z[p]
    xcm_sys = xcm_sys / total_particles
    ycm_sys = ycm_sys / total_particles
    zcm_sys = zcm_sys / total_particles
    if(abs(xcm_sys-xcm_sys_zero) > 0.01 or \
       abs(ycm_sys-ycm_sys_zero) > 0.01 or \
       abs(zcm_sys-zcm_sys_zero) > 0.01):
      sys.stdout.write('r_COM(t=%d steps): %8.2f %8.2f %8.2f\n' % \
      (timestep,xcm_sys-xcm_sys_zero,ycm_sys-ycm_sys_zero,zcm_sys-zcm_sys_zero))

    # correct positions
    if 0:
      for i in range(M_rings + M_linear):
	for j in range(monomers):
	  p = i * monomers + j
	  x[p] = x[p] - (xcm_sys - xcm_sys_zero)
	  y[p] = y[p] - (ycm_sys - ycm_sys_zero)
	  z[p] = z[p] - (zcm_sys - zcm_sys_zero)

    # initialize ring values
    g1_rings = 0.0
    g2_rings = 0.0
    g3_rings = 0.0

    # loop over rings
    for i in range(M_rings):

      # compute center-of-mass
      xcm = 0.0; ycm = 0.0; zcm = 0.0
      for j in range(monomers):
	p = i * monomers + j
	xcm = xcm + x[p]
	ycm = ycm + y[p]
	zcm = zcm + z[p]
	g1_rings = g1_rings + (x[p] - xi[p])**2 + (y[p] - yi[p])**2 + (z[p] - zi[p])**2
      xcm = xcm / monomers
      ycm = ycm / monomers
      zcm = zcm / monomers

      # compute g2
      for j in range(monomers):
	p = i * monomers + j
	g2_rings = g2_rings + (x[p] - xi[p] - xcm + com_zero[i][0])**2 + \
		              (y[p] - yi[p] - ycm + com_zero[i][1])**2 + \
		              (z[p] - zi[p] - zcm + com_zero[i][2])**2

      # compute g3
      g3_rings = g3_rings + (xcm - com_zero[i][0])**2 \
                          + (ycm - com_zero[i][1])**2 \
                          + (zcm - com_zero[i][2])**2

    # initialize linear values
    g1_linear = 0.0
    g1_linear_inner = 0.0
    g2_linear = 0.0
    g3_linear = 0.0

    # loop over linear chains
    for i in range(M_linear):

      # compute center-of-mass
      xcm = 0.0; ycm = 0.0; zcm = 0.0
      for j in range(monomers):
        p = M_rings * monomers + i * monomers + j
        xcm = xcm + x[p]
        ycm = ycm + y[p]
        zcm = zcm + z[p]
        g1_linear = g1_linear + (x[p] - xi[p])**2 + (y[p] - yi[p])**2 + (z[p] - zi[p])**2
      xcm = xcm / monomers
      ycm = ycm / monomers
      zcm = zcm / monomers

      jstart = monomers_half - inner_monomers / 2
      jend = jstart + inner_monomers
      mnrs = range(jstart, jend)
      if(len(mnrs) != inner_monomers): print "inner monomers"; sys.exit(1)
      for j in mnrs:
        p = M_rings * monomers + i * monomers + j
        g1_linear_inner = g1_linear_inner + (x[p] - xi[p])**2 + (y[p] - yi[p])**2 + (z[p] - zi[p])**2

      # compute g2
      for j in range(monomers):
        p = M_rings * monomers + i * monomers + j
        g2_linear = g2_linear + (x[p] - xi[p] - xcm + com_zero[M_rings + i][0])**2 + \
                                (y[p] - yi[p] - ycm + com_zero[M_rings + i][1])**2 + \
                                (z[p] - zi[p] - zcm + com_zero[M_rings + i][2])**2

      # compute g3
      g3_linear = g3_linear + (xcm - com_zero[M_rings + i][0])**2 \
                            + (ycm - com_zero[M_rings + i][1])**2 \
                            + (zcm - com_zero[M_rings + i][2])**2

    # store time difference
    dt = timestep - time_zero

    if(not cntr.has_key(dt)):
      cntr[dt] = 0
      g1ave_rings[dt] = 0.0
      g2ave_rings[dt] = 0.0
      g3ave_rings[dt] = 0.0
      g1ave_linear[dt] = 0.0
      g1ave_linear_inner[dt] = 0.0
      g2ave_linear[dt] = 0.0
      g3ave_linear[dt] = 0.0
 
    cntr[dt] = cntr[dt] + 1
    if(M_rings != 0):
      g1ave_rings[dt] = g1ave_rings[dt] + g1_rings / (M_rings * monomers)
      g2ave_rings[dt] = g2ave_rings[dt] + g2_rings / (M_rings * monomers)
      g3ave_rings[dt] = g3ave_rings[dt] + g3_rings / M_rings
    if(M_linear != 0):
      g1ave_linear[dt] = g1ave_linear[dt] + g1_linear / (M_linear * monomers)
      g1ave_linear_inner[dt] = g1ave_linear_inner[dt] + g1_linear_inner / (M_linear * inner_monomers)
      g2ave_linear[dt] = g2ave_linear[dt] + g2_linear / (M_linear * monomers)
      g3ave_linear[dt] = g3ave_linear[dt] + g3_linear / M_linear

# normalize data
rings = []
linear = []
for dt in cntr.keys():
  rings.append((dt, g1ave_rings[dt] / cntr[dt], g2ave_rings[dt] / cntr[dt], g3ave_rings[dt] / cntr[dt]))
  linear.append((dt, g1ave_linear[dt] / cntr[dt], g2ave_linear[dt] / cntr[dt], g3ave_linear[dt] / cntr[dt], g1ave_linear_inner[dt] / cntr[dt]))

# sort by time
rings.sort(lambda u, v: cmp(u[0], v[0]))
linear.sort(lambda u, v: cmp(u[0], v[0]))

# outfile 
MringsMlinear = ''
if(M_linear != 0): MringsMlinear = str(M_rings) + '_' + str(M_linear)

if(M_rings != 0):
  # write rings data to file
  outfile = 'g123_' + MringsMlinear + 'N' + str(monomers) + case + '_rings.dat'
  fout = open(outfile, 'w')
  fout.write('# ' + time.asctime() + '\n')
  fout.write('# ' + os.getcwd() + '\n')
  fout.write('# M_rings  M_linear  monomers\n')
  fout.write('# ' + str(M_rings) + '  ' + str(M_linear) + '  ' + str(monomers) + '\n')
  fout.write('# time_step  g1  g2  g3 (sigma**2)\n')
  for item in rings:
    fout.write('%10d %10.3f %10.3f %10.3f\n' % (item[0], item[1], item[2], item[3]))
  fout.close()
  sys.stdout.write(outfile + ' has been written to disk.\n')

if(M_linear != 0):
  # write linear chain data to file
  outfile = 'g123_' + MringsMlinear + 'N' + str(monomers) + case + '_linear.dat'
  fout = open(outfile, 'w')
  fout.write('# ' + time.asctime() + '\n')
  fout.write('# ' + os.getcwd() + '\n')
  fout.write('# M_rings  M_linear  monomers\n')
  fout.write('# ' + str(M_rings) + '  ' + str(M_linear) + '  ' + str(monomers) + '\n')
  fout.write('# time_step  g1  g2  g3  g1_in (sigma**2)\n')
  for item in linear:
    fout.write('%10d %10.3f %10.3f %10.3f %10.3f\n' % (item[0], item[1], item[2], item[3], item[4]))
  fout.close()
  sys.stdout.write(outfile + ' has been written to disk.\n')

"""This script computes Rg squared and Re squared for a system
   of rings and/or linear chains. The average values and
   standard deviations are reported. Coordinates are assumed
   to be unwrapped."""

# use sys for printing compatibility between Python 2 & 3
import sys
import time
import os
import math
import glob
import pypolymer

case = ''

# generation of file names: (1) glob, (2) loop, (3) manual or (4) log scale
blockfiles_method = 1

if(blockfiles_method == 1):

  prefix = '../../'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]
  if 1:
    maxfiles = 25
    numfiles = len(blockfiles)
    if(numfiles > maxfiles): blockfiles = blockfiles[::int(numfiles / float(maxfiles))]

elif(blockfiles_method == 2):

  prefix = '../../'
  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append(prefix + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['../../']

elif(blockfiles_method == 4):

  # prefix
  prefix = '../../'

  # number of files per decade
  n = 5

  # exponent of first and last decade (time steps)
  start_exp = 5
  end_exp = 10
 
  # base
  b = math.e

  number_of_decades = end_exp - start_exp
  start = math.log(10**start_exp)
  end = math.log(10**end_exp)
  step = (end - start) / (n * number_of_decades - 1)
  uniform_time = []
  for i in range(n * number_of_decades):
    uniform_time.append(b**(start + step * i))
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles = [int(u[u.rindex('.') + 1:]) for u in existingfiles]
  existingfiles.sort()
  uniform_time = [pypolymer.find_closest_existing_file(existingfiles, u) for u in uniform_time]
  uniform_time = pypolymer.unique(uniform_time)
  blockfiles = [prefix + str(u) for u in uniform_time]

else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit(1)

if(blockfiles == []):
  sys.stderr.write('Blockfiles is empty. Exiting ...\n')
  sys.exit(1)

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                      Rg2 Re2                    \n')
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

# initialize the lists (time, Rgsq, Resq)
time_rsq_rings = []
time_rsq_linear = []

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
    sys.exit(1)

  sys.stdout.write('Working on ' + file + '\n')

  if(M_rings != 0):
    # loop over rings
    radsq = 0.0
    spansq = 0.0
    for i in range(M_rings):

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
      for j in range(monomers):
	p = i * monomers + j
	radsq = radsq + (x[p] - xcm)**2 + (y[p] - ycm)**2 + (z[p] - zcm)**2

      # compute spanning distances squared
      for j in range(monomers_half):
	p1 = i * monomers + j
	p2 = p1 + monomers_half
	spansq = spansq + (x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1] - z[p2])**2

    # compute radius of gyration squared
    RGsq = radsq / (M_rings * monomers)

    # compute diameter squared
    REsq = spansq / (M_rings * monomers_half)

    # append to list
    time_rsq_rings.append((timestep, RGsq, REsq))

  if(M_linear != 0):
    # loop over linear chains
    radsq = 0.0
    spansq = 0.0
    for i in range(M_linear):

      # compute center-of-mass
      xcm = 0.0; ycm = 0.0; zcm = 0.0
      for j in range(monomers):
	p = M_rings * monomers + i * monomers + j
	xcm = xcm + x[p]
	ycm = ycm + y[p]
	zcm = zcm + z[p]
      xcm = xcm / monomers
      ycm = ycm / monomers
      zcm = zcm / monomers

      # sum the squares
      for j in range(monomers):
	p = M_rings * monomers + i * monomers + j
	radsq = radsq + (x[p] - xcm)**2 + (y[p] - ycm)**2 + (z[p] - zcm)**2

      # compute end-to-end distance squared
      p1 = M_rings * monomers + i * monomers
      p2 = p1 + monomers - 1
      spansq = spansq + (x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1] - z[p2])**2

    # compute radius of gyration squared
    RGsq = radsq / (M_linear * monomers)

    # compute end-to-end distance squared
    REsq = spansq / M_linear

    # append to list
    time_rsq_linear.append((timestep, RGsq, REsq))

# outfile 
MringsMlinear = ''
if(M_linear != 0): MringsMlinear = str(M_rings) + '_' + str(M_linear)

if(M_rings != 0):

  # compute the averages
  ave_RGsq = 0.0
  ave_REsq = 0.0
  for i in range(len(time_rsq_rings)):
    ave_RGsq = ave_RGsq + time_rsq_rings[i][1]
    ave_REsq = ave_REsq + time_rsq_rings[i][2]
  ave_RGsq = ave_RGsq / len(time_rsq_rings)
  ave_REsq = ave_REsq / len(time_rsq_rings)

  # compute the standard deviations
  Gsq = 0
  Esq = 0
  for i in range(len(time_rsq_rings)):
    Gsq = Gsq + (time_rsq_rings[i][1] - ave_RGsq)**2
    Esq = Esq + (time_rsq_rings[i][2] - ave_REsq)**2
  RGsq_error = math.sqrt(Gsq / len(time_rsq_rings))
  REsq_error = math.sqrt(Esq / len(time_rsq_rings))

  # write data to file
  outfile = 'RG2RE2_' + MringsMlinear + 'N' + str(monomers) + case + '_rings.dat'
  fout = open(outfile, 'w')
  fout.write('# ' + time.asctime() + '\n')
  fout.write('# ' + os.getcwd() + '\n')
  fout.write('# M_rings M_linear monomers\n')
  fout.write('# ' + str(M_rings) + ' ' + str(M_linear) + ' ' + str(monomers) + '\n')
  fout.write('# RGsq is the radius of gyration squared\n')
  fout.write('# REsq is the diameter squared\n#\n')
  fout.write('# time in tau\n')
  fout.write('# ave_RGsq, RGsq_error, ave_REsq and REsq_error in sigma**2\n')
  fout.write('# RGsq_error and REsq_error correspond to one standard deviation.\n#\n')
  fout.write('# ave_RGsq RGsq_error ave_REsq REsq_error\n')
  fout.write('# %8.2f %8.2f %8.2f %8.2f\n' % (ave_RGsq, RGsq_error, ave_REsq, REsq_error))
  fout.write('# time (tau)  RGsq (sigma**2)  REsq (sigma**2)\n')
  for rec in time_rsq_rings:
    fout.write('%10d %8.2f %8.2f\n' % (rec[0], rec[1], rec[2]))
  fout.close()
  sys.stdout.write(outfile + ' has been written to disk.\n')

if(M_linear != 0):

  # compute the averages
  ave_RGsq = 0.0
  ave_REsq = 0.0
  for i in range(len(time_rsq_linear)):
    ave_RGsq = ave_RGsq + time_rsq_linear[i][1]
    ave_REsq = ave_REsq + time_rsq_linear[i][2]
  ave_RGsq = ave_RGsq / len(time_rsq_linear)
  ave_REsq = ave_REsq / len(time_rsq_linear)

  # compute the standard deviations
  Gsq = 0
  Esq = 0
  for i in range(len(time_rsq_linear)):
    Gsq = Gsq + (time_rsq_linear[i][1] - ave_RGsq)**2
    Esq = Esq + (time_rsq_linear[i][2] - ave_REsq)**2
  RGsq_error = math.sqrt(Gsq / len(time_rsq_linear))
  REsq_error = math.sqrt(Esq / len(time_rsq_linear))

  # write data to file
  outfile = 'RG2RE2_' + MringsMlinear + 'N' + str(monomers) + case + '_linear.dat'
  fout = open(outfile, 'w')
  fout.write('# ' + time.asctime() + '\n')
  fout.write('# ' + os.getcwd() + '\n')
  fout.write('# M_rings M_linear monomers\n')
  fout.write('# ' + str(M_rings) + ' ' + str(M_linear) + ' ' + str(monomers) + '\n')
  fout.write('# RGsq is the radius of gyration squared\n')
  fout.write('# REsq is the end-to-end distance squared\n#\n')
  fout.write('# time in tau\n')
  fout.write('# ave_RGsq, RGsq_error, ave_REsq and REsq_error in sigma**2\n')
  fout.write('# RGsq_error and REsq_error correspond to one standard deviation.\n#\n')
  fout.write('# ave_RGsq RGsq_error ave_REsq REsq_error\n')
  fout.write('# %8.2f %8.2f %8.2f %8.2f\n' % (ave_RGsq, RGsq_error, ave_REsq, REsq_error))
  fout.write('# time (tau)  RGsq (sigma**2)  REsq (sigma**2)\n')
  for rec in time_rsq_linear:
    fout.write('%10d %8.2f %8.2f\n' % (rec[0], rec[1], rec[2]))
  fout.close()
  sys.stdout.write(outfile + ' has been written to disk.\n')

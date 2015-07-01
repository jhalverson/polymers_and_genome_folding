"""This script computes the center-of-mass position and
   velocity of the system."""

# use sys for printing compatibility between Python 2 & 3
import sys
import os
import math

# all systems have two hundred chains
chains = 200

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

elif(blockfiles_method == 2):

  istr = 370000000
  iend = 380000000
  incr = 5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M800a.013160000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

# initialize the list (time, position, velocity)
rv = []

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

  # store the coordinates and velocities
  x = []; y = []; z = []
  vx = []; vy = []; vz = []
  for i in range(total_particles):
    sline = f.readline().split()
    x.append(float(sline[1]))
    y.append(float(sline[2]))
    z.append(float(sline[3]))
    vx.append(float(sline[5]))
    vy.append(float(sline[6]))
    vz.append(float(sline[7]))

  # close the file
  f.close()

  # initialize COM position of the system
  xcm_sys = 0.0
  ycm_sys = 0.0
  zcm_sys = 0.0

  # compute COM position
  for i in range(chains):
    for j in range(monomers):
      p = i * monomers + j
      xcm_sys = xcm_sys + x[p]
      ycm_sys = ycm_sys + y[p]
      zcm_sys = zcm_sys + z[p]
  xcm_sys = xcm_sys / total_particles
  ycm_sys = ycm_sys / total_particles
  zcm_sys = zcm_sys / total_particles

  # initialize COM velocity of the system
  vxcm_sys = 0.0
  vycm_sys = 0.0
  vzcm_sys = 0.0

  # compute COM position
  for i in range(chains):
    for j in range(monomers):
      p = i * monomers + j
      vxcm_sys = vxcm_sys + vx[p]
      vycm_sys = vycm_sys + vy[p]
      vzcm_sys = vzcm_sys + vz[p]
  vxcm_sys = vxcm_sys / total_particles
  vycm_sys = vycm_sys / total_particles
  vzcm_sys = vzcm_sys / total_particles

  # compute magnitude of COM velocity
  vcm_magn = math.sqrt(vxcm_sys**2 + vycm_sys**2 + vzcm_sys**2)

  # store values
  rv.append((timeval, xcm_sys, ycm_sys, zcm_sys, vxcm_sys, vycm_sys, vzcm_sys, vcm_magn))

# write out data
outfile = 'COM' + str(monomers) + '.dat'
fdat = open(outfile, 'w')
fdat.write('# time (tau) xcm ycm zcm (sigma) vxcm vycm vzcm vmagn (sigma/tau)\n')
for rec in rv:
  t = rec[0]
  x = rec[1]
  y = rec[2]
  z = rec[3]
  u = rec[4]
  v = rec[5]
  w = rec[6]
  m = rec[7]
  fdat.write('%10d %.2f %.2f %.2f %g %g %g %g\n' % (t, x, y, z, u, v, w, m))
fdat.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

# write a Gnuplot input file
pltfile = 'COM' + str(monomers) + '.plt'
f = open(pltfile, 'w')
f.write('set xlabel "TIME (TAU)" \n')
f.write('set ylabel "VELOCITY (SIGMA/TAU)" \n')
f.write('plot \\\n')
f.write('"' + outfile + '" using 1:5 with lines title "X-VELOCITY", \\\n')
f.write('"' + outfile + '" using 1:6 with lines title "Y-VELOCITY", \\\n')
f.write('"' + outfile + '" using 1:7 with lines title "Z-VELOCITY"')
f.close()
sys.stdout.write(pltfile + ' has been written to disk (Gnuplot input file).\n')

if(gnuplot_load == 'yes'):
  cmd = 'gnuplot -persist ' + pltfile
  os.system(cmd)

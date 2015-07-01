"""This script tracks the distance between particles and their
   individual motion.

   1. find the 200 COM positions in the central cell
   2. translate these to the 26 surrounding cells
   3. this produces 5400 chains (we do not need to consider others)
   4. store all the COM positions
   5. read in a second file with a later time value
   6. the central chains have new positions of r_new = r_0 + r(t) - r(0)
   7. these can be translated to the neighboring cells
   8. distances can be computed
"""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os

# all systems have two hundred chains
chains = 200
chains26 = 26 * chains

# density
rho = 0.85

# intermediate output ('yes' or 'no')
verbose = 'yes'

# load data into Gnuplot ('yes' or 'no')
gnuplot_load = 'no'

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 2

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('*000000')

  # sort blockfiles list
  blockfiles.sort(lambda u, v: cmp(int(u[u.index('.') + 1:]), int(v[v.index('.') + 1:])))

elif(blockfiles_method == 2):

  istr = 140000000
  iend = 440000000
  incr =   5000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RGEQ_N200M400.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RGEQ_N200M400.0130700000', 'RGEQ_N200M400.0200000000', 'RGEQ_N200M400.0439880000']
 
else:

  sys.stderr.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')

# wrapped COM positions at time origin
com_wrap_zero = []

# unwrapped COM positions at time origin
com_unwrap_zero = []

# unwrapped COM positions at time origin
com_wrap_trans_zero = []

# list of sep3's
master = []

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
  monomers = total_particles / chains

  # safety checks
  if(file == blockfiles[0]):
    monomers1st = monomers
    sys.stdout.write('chains = ' + str(chains) + ', monomers = ' + str(monomers) + '\n\n')
  if(monomers != monomers1st):
    sys.stderr.write('ERROR: monomers is not equal to monomers1st\n')
    sys.exit(1)
  if(verbose == 'yes'):
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
    r = f.readline().split()
    x.append(float(r[1]))
    y.append(float(r[2]))
    z.append(float(r[3]))

  # close the file
  f.close()

  # use first file to get initial positions
  if(file == blockfiles[0]):
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

      # store unwrap COM positions at t = 0
      com_unwrap_zero.append((xcm, ycm, zcm))

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
      com_wrap_zero.append((xcm, ycm, zcm))

    # translate and store the wrapped COM positions
    for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          if(i**2 + j**2 + k**2 != 0):
            for c in range(chains):
              xcm = com_wrap_zero[c][0] + i * L
              ycm = com_wrap_zero[c][1] + j * L
              zcm = com_wrap_zero[c][2] + k * L
              com_wrap_trans_zero.append((xcm, ycm, zcm))
 
  # if not first file
  else:
    com_unwrap_t = []
    com_central_t = []
    com_trans_t = []

    # compute COM at time t
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

      # store center-of-mass position at time t
      com_unwrap_t.append((xcm, ycm, zcm))

    # get COM positions of central chains at time t
    for i in xrange(chains):
      xcm = com_wrap_zero[i][0] + (com_unwrap_t[i][0] - com_unwrap_zero[i][0])
      ycm = com_wrap_zero[i][1] + (com_unwrap_t[i][1] - com_unwrap_zero[i][1])
      zcm = com_wrap_zero[i][2] + (com_unwrap_t[i][2] - com_unwrap_zero[i][2])
      com_central_t.append((xcm, ycm, zcm))

    # translate and store the wrapped COM positions
    m = 0
    for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          if(i**2 + j**2 + k**2 != 0):
            for c in xrange(chains):
              xcm = com_wrap_trans_zero[m][0] + (com_unwrap_t[c][0] - com_unwrap_zero[c][0])
              ycm = com_wrap_trans_zero[m][1] + (com_unwrap_t[c][1] - com_unwrap_zero[c][1])
              zcm = com_wrap_trans_zero[m][2] + (com_unwrap_t[c][2] - com_unwrap_zero[c][2])
              com_trans_t.append((xcm, ycm, zcm))
              m = m + 1

    # find minimum distance between central and all others
    for i in xrange(chains):
      sep = []
      sep3 = []
      for j in xrange(chains):
        if(i != j):
          xij = com_central_t[i][0] - com_central_t[j][0]
          yij = com_central_t[i][1] - com_central_t[j][1]
          zij = com_central_t[i][2] - com_central_t[j][2]
          rijsq = xij**2 + yij**2 + zij**2
          sep.append((i, j, rijsq))

      for j in xrange(chains26):
        xij = com_central_t[i][0] - com_trans_t[j][0]
        yij = com_central_t[i][1] - com_trans_t[j][1]
        zij = com_central_t[i][2] - com_trans_t[j][2]
        rijsq = xij**2 + yij**2 + zij**2
        sep.append((i, j + chains, rijsq))

      # sort list
      sep.sort(lambda u, v: cmp(u[2], v[2]))
      
      # extract three closest
      sep3 = [sep[0], sep[1], sep[2]]
      #sep3.sort(lambda u, v: cmp(u[1], v[1]))
      master.append(sep3)

      fname = str(i) + '.ngh'
      f = open(fname, 'a')
      f.write('%10d %10g %10g %10g\n' % (timeval, sep3[0][1], sep3[1][1], sep3[2][1]))
      f.close()

# perform analysis on sep3's
master_pairs = []
rows = len(blockfiles) - 1
for i in range(chains):
  ct = []
  for a in xrange(5399):
    ct.append(0)
  for j in range(rows):
    s = i + chains * j
    sep3 = master[s]
    for k in [0, 1, 2]:
      n = sep3[k][1]
      ct[n] = ct[n] + 1

  pair = []
  for a in xrange(5399):
    pair.append((a, ct[a]))
  # reverse sort
  pair.sort(lambda u, v: cmp(u[1], v[1]), reverse=True)
  master_pairs.append((i, pair[0], pair[1], pair[2]))
  print i, pair[0], pair[1], pair[2]

# reverse sort
master_pairs.sort(lambda u, v: cmp(u[1][1], v[1][1]), reverse=True)
for i in range(chains):
  print master_pairs[i]

sys.exit()

# write Gnuplot script
fname = 'correlated_diffusion.plt'
fgnu = open(fname, 'w')
fgnu.write('set terminal png medium\n')
fgnu.write('set xlabel "TIME (TAU)"\n')
fgnu.write('set ylabel "NEIGHBOR ID"\n')
fgnu.write('set yrange [0:5399]\n')
for i in range(chains):
  fgnu.write('set title "THREE CLOSEST NEIGHBORS FOR CHAIN ' + str(i) + '"\n')
  fgnu.write('set output "' + str(i) + '.png"\n')
  fin = '"' + str(i) + '.ngh"'
  fgnu.write('plot ' + fin + ' using 1:2, ' + fin + ' using 1:3, ' + fin + ' using 1:4\n\n')
fgnu.close()

# call Gnuplot
err = os.system('gnuplot correlated_diffusion.plt')

# create HTML page
fname = 'correlated_diffusion.html'
fhtml = open(fname, 'w')
fhtml.write('<html>\n')
fhtml.write('<head><title>' + str(monomers) + 'a</title></head>\n')
fhtml.write('<body>\n')
for i in range(chains):
  img = '"' + str(i) + '.png"'
  fhtml.write('<img src=' + img + ' border="1"><p>\n')
fhtml.write('</body>\n')
fhtml.write('</html>')
fhtml.close()

"""This script computes two orthogonal spanning vectors of each chain
   by averaging the positions of n monomers in forming the
   head and tail point of each vector. The cross product of
   these vectors is also computed. Integration time step is assumed
   to be 0.01 tau."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import glob
import pypolymer

# all systems have two hundred chains
chains = 200

# how many monomers are averaged to form a point
n = 5

# simulation run a, b or c (this is used to select the shift)
run = 'c'

# generation of file names: (1) glob, (2) loop, (3) manual
blockfiles_method = 2

if(blockfiles_method == 1):

  blockfiles = glob.glob('*000000')
  blockfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), blockfiles)
  blockfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 500000
  iend = 2000000000
  incr = 500000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('../stripped/' + str(i) + '.pos')

elif(blockfiles_method == 3):

  blockfiles = ['']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n\n')

# initialize list
v = []

# estimate of ReSq and the angle between a and b
spanSq = 0.0
angle = 0.0

# loop over block files
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
    shift = {'a':0, 'b':monomers/12, 'c':monomers/6}
    shift = shift[run]
  if(total_particles1st != total_particles):
    sys.stderr.write('Change in total_particles (' + str(total_particles) + ').\n')
    sys.exit(1)

  sys.stdout.write('Working on ' + file + '\n')

  # initialize sublist
  v_sub = []

  # loop over chains
  for i in range(chains):

    x1head = 0.0
    y1head = 0.0
    z1head = 0.0

    x1tail = 0.0
    y1tail = 0.0
    z1tail = 0.0

    for j in range(n):

      p = i * monomers + j + shift
      x1head = x1head + x[p]
      y1head = y1head + y[p]
      z1head = z1head + z[p]

      p = i * monomers + j + monomers / 2 + shift
      x1tail = x1tail + x[p]
      y1tail = y1tail + y[p]
      z1tail = z1tail + z[p]

    x1head = x1head / n
    y1head = y1head / n
    z1head = z1head / n

    x1tail = x1tail / n
    y1tail = y1tail / n
    z1tail = z1tail / n

    # store vector a
    a = (x1head - x1tail, y1head - y1tail, z1head - z1tail)

    x2head = 0.0
    y2head = 0.0
    z2head = 0.0

    x2tail = 0.0
    y2tail = 0.0
    z2tail = 0.0

    for j in range(n):

      p = i * monomers + j + monomers / 4 + shift
      x2head = x2head + x[p]
      y2head = y2head + y[p]
      z2head = z2head + z[p]

      p = i * monomers + j + 3 * monomers / 4 + shift
      x2tail = x2tail + x[p]
      y2tail = y2tail + y[p]
      z2tail = z2tail + z[p]

    x2head = x2head / n
    y2head = y2head / n
    z2head = z2head / n

    x2tail = x2tail / n
    y2tail = y2tail / n
    z2tail = z2tail / n

    # store vector b
    b = (x2head - x2tail, y2head - y2tail, z2head - z2tail)

    # compute spanning distance squared and angle
    magn_a_sq = a[0]**2 + a[1]**2 + a[2]**2
    magn_b_sq = b[0]**2 + b[1]**2 + b[2]**2
    spanSq = spanSq + magn_a_sq
    spanSq = spanSq + magn_b_sq
    magn_a = magn_a_sq**0.5
    magn_b = magn_b_sq**0.5
    angle = angle + math.acos((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (magn_a * magn_b)) * (180.0 / math.pi)

    # compute cross product
    c = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

    v_sub.append([timestep, a, b, c])

  # append list with (timestep, a, b, c) for each chain
  v.append(v_sub)

# write out the data
outfile = 'ab' + str(monomers) + run + '_n' + str(n) + '.dat'
fmt = '%d %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n'
fout = open(outfile, 'w')
fout.write('# time (tau)  chain_id  a_x  a_y  a_z  b_x  b_y  b_z (sigma)  c_x  c_y  c_z (sigma**2)\n')
for k in range(len(blockfiles)):
  for i in range(chains):
    a = v[k][i][1]
    b = v[k][i][2]
    c = v[k][i][3]
    outline = fmt % (v[k][i][0] / 100, i, a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2])
    fout.write(outline)
fout.close()

sys.stdout.write(outfile + ' has been written to disk.\n')
sys.stdout.write('Estimated ReSq = %.2f\n' % (spanSq / (2 * chains * len(blockfiles))))
sys.stdout.write('Average angle = %.1f\n' % (angle / (chains * len(blockfiles))))

"""This Python script computes the zero-shear viscosity by
   numerically integrating the time-dependent shear modulus G(t)."""

import sys
import os
import time

# input file
fin = 'grest_acf_M250N200_all.dat'

# bin width (tau)
dt = 100 * 0.012

# read and store file contents
f = open(fin)
lines = f.readlines()
f.close()

# load data
gt = []
for i in range(len(lines)):
  line = lines[i]
  s = line.split()
  gt.append(float(s[1]))

sum_gtdt = 0.5 * dt * gt[0]
viscosity = [(0.0, sum_gtdt)]
for i in xrange(1, len(gt)):
  sum_gtdt = sum_gtdt + dt * gt[i]
  if(i % 1000 == 0): viscosity.append((i * dt, sum_gtdt))

outfile = fin + '.viscosity'
f = open(outfile, 'w')
f.write('# ' + time.asctime() + '\n')
f.write('# ' + os.getcwd() + '\n')
f.write('# t (tau)  viscosity (epsilon tau / sigma**3)\n')
for i in xrange(len(viscosity)):
  f.write('%.1f %.6e\n' % (viscosity[i][0], viscosity[i][1]))
f.close()
sys.stdout.write(outfile + ' has been written.\n')

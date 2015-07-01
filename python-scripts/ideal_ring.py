"""This script generates n equally-spaced points on a circle.
   Randomness may be applied. The coordinates are written
   to a PDB file."""

import math
import random

# number of points
n = 100

# radius of ring
radius = 17.0

# PDB format string
pdb = '%6s%5d%5s%4s%6d%12.3f%8.3f%8.3f\n'

# scale factor for coordinates
sf = 2.0

# scale factor for random numbers
rsf = 0.25

# initialize lists
x = []
y = []
z = []

# function which returns scaled random numbers
def rnd():
  return rsf * random.uniform(-1.0, 1.0)

# compute coordinates
for i in range(n):
  step = 2.0 * math.pi / n
  theta = 0.0 + i * step
  x.append(radius * math.cos(theta) + rnd())
  y.append(radius * math.sin(theta) + rnd())
  z.append(rnd())

# write data to PDB file
fname = 'ideal_ring' + str(n) + '.pdb'
fpdb = open(fname, 'w')
for i in range(n):
  fpdb.write(pdb % ('HETATM', i, 'C', 'ATM', i, sf * x[i], sf * y[i], sf * z[i]))
print '\n' + fname + ' has been written to disk.\n'

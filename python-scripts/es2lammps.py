"""This script converts an ESPResSo block file for a polymer 
   melt of rings to a LAMMPS input file."""

import sys

# number of chains
chains = 200

# file to convert
f = open('RGEQ_N200M400.0548200000')

# find line where box size is found and store it
line = ''
while not '{box_l ' in line:
  line = f.readline()
box = line.replace('}','').split()
xlen, ylen, zlen = float(box[1]), float(box[2]), float(box[3])
sys.stdout.write('(Lx, Ly, Lz) = (%.4f, %.4f, %.4f)\n' % (xlen, ylen, zlen))

# find line where total particles is found and store it
line = ''
while not '{n_part ' in line:
  line = f.readline()
np = line.replace('}','').split()
total_particles = int(np[1])
monomers = total_particles / chains
sys.stdout.write('(monomers, chains, total particles) = (%d, %d, %d)\n' % \
                (monomers, chains, total_particles))

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

# wrap chains by COM
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

  # shift chains such that the central box is between [0, L]
  ix = int(xcm / xlen)
  iy = int(ycm / ylen)
  iz = int(zcm / zlen)

  # correct in case of negative positions
  if(xcm < 0.0): ix = ix - 1
  if(ycm < 0.0): iy = iy - 1
  if(zcm < 0.0): iz = iz - 1

  # shift coordinates so that center-of-mass in in central box
  for j in range(monomers):
    p = i * monomers + j
    x[p] = x[p] - ix * xlen
    y[p] = y[p] - iy * ylen
    z[p] = z[p] - iz * zlen

# find line where topology begins
line = ''
while not '{bonds' in line:
  line = f.readline()

# store the bond information
bonds = []; angles = []
for i in range(total_particles):
  line = f.readline()
  line = line.replace('{','')
  bnd =  line.replace('}','').split()
  bonds.append((int(bnd[0]), int(bnd[2])))
  if(len(bnd) == 6):
    angles.append((int(bnd[0]), int(bnd[4]), int(bnd[5])))

f.close()

# write out LAMMPS input file
top_definition = """LAMMPS

80000 atoms
80000 bonds
80000 angles
    0 dihedrals
    0 impropers

    1 atom types
    1 bond types
    1 angle types
    0 dihedral types
    0 improper types
"""

flnm = '400a.rings'
F = open(flnm, 'w')
F.write(top_definition)
F.write('\n')
F.write('0.0000 %.4f xlo xhi\n' % (xlen))
F.write('0.0000 %.4f ylo yhi\n' % (ylen))
F.write('0.0000 %.4f zlo zhi\n' % (zlen))
F.write('\n')
F.write('Masses\n')
F.write('\n')
F.write('1 1.0\n')
F.write('\n')

F.write('Atoms\n\n')
for i in range(total_particles):
  molecule = int(i / monomers) + 1
  outline = '%6d %6d %6d %12.8f %12.8f %12.8f\n' % (i + 1, molecule, 1, x[i], y[i], z[i])
  #outline = '%6d %6d %12.8f %12.8f %12.8f\n' % (i + 1, 1, x[i], y[i], z[i])
  F.write(outline)

F.write('\nBonds\n\n')
for i in range(len(bonds)):
  outline = '%6d %6d %6d %6d\n' % (i + 1, 1,  bonds[i][1] + 1, bonds[i][0] + 1)
  F.write(outline)

F.write('\nAngles\n\n')
for i in range(len(angles)):
  outline = '%6d %6d %6d %6d %6d\n' % (i + 1, 1, angles[i][1] + 1, angles[i][0] + 1, angles[i][2] + 1)
  F.write(outline)

F.close()
sys.stdout.write('LAMMPS input file has been written to disk.\n')

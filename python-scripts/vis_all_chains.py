"""This script writes all chains to a PDB file. The
   coordinates are wrapped to those of the central
   simulation cell. A PDB file of the box may also
   be written out."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import os

# all systems have two hundred chains
chains = 200

# density
rho = 0.85

# display intermediate output ('yes' or 'no')
verbose = 'yes'

# load output into PyMol ('yes' or 'no')
pymol_load = 'no'

# PDB format string
pdb = '%6s%5d%5s%4s%6d%12.3f%8.3f%8.3f\n'

# scale factor for coordinates
sf = 2.0

# write box file ('yes' or 'no')
write_box = 'yes'

# generation of file names: (1) glob, (2) loop, (3) manual
blockfiles_method = 3

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('*000000')

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

  blockfiles = ['RG_N200M1600.0361200000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

sys.stdout.write('\n\nTotal number of block files: ' + str(len(blockfiles)) + '\n')

# loop over block files
for file in blockfiles:

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

  # loop over chains to wrap coordinates
  for i in range(chains):

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

    # shift chain such that the central box is between [0, L]
    ix = int(xcm / L)
    iy = int(ycm / L)
    iz = int(zcm / L)
 
    # correct in case of negative positions
    if(xcm < 0.0): ix = ix - 1
    if(ycm < 0.0): iy = iy - 1
    if(zcm < 0.0): iz = iz - 1
  
    # perform wrap
    for j in range(monomers):
      p = i * monomers + j
      x[p] = x[p] - ix * L
      y[p] = y[p] - iy * L
      z[p] = z[p] - iz * L

  # determine files and limits
  pdb_max_particles = 99999
  max_chains_per_file = int(pdb_max_particles / monomers)
  number_of_full_files = int(chains / max_chains_per_file)
  chains_in_final_file = chains - number_of_full_files * max_chains_per_file

  # determine if final file needs to be written
  final_file = 1
  if(chains_in_final_file == 0): final_file = 0

  # loop over full files
  for k in range(number_of_full_files):
    fname = 'M' + str(monomers) + 'T' + str(timeval) + 'F' + str(k) + '.pdb'
    outfile = open(fname, 'w')
    for i in range(max_chains_per_file):
      for j in range(monomers):
        chain_id = k * max_chains_per_file + i
        p = chain_id * monomers + j
        q = i * monomers + j
        outfile.write(pdb % ('HETATM', q, 'C', 'ATM', chain_id, sf * x[p], sf * y[p], sf * z[p]))
    outfile.close()
    sys.stdout.write('  ' + fname + ' written to disk.\n')

  # write out remaining chains (if any) to final file
  if(final_file == 1):
    fname = 'M' + str(monomers) + 'T' + str(timeval) + 'F' + str(number_of_full_files) + '.pdb'
    outfile = open(fname, 'w')
    for i in range(chains_in_final_file):
      for j in range(monomers):
        chain_id = number_of_full_files * max_chains_per_file + i
        p = chain_id * monomers + j
        q = i * monomers + j
        outfile.write(pdb % ('HETATM', q, 'C', 'ATM', chain_id, sf * x[p], sf * y[p], sf * z[p]))
    outfile.close()
    sys.stdout.write('  ' + fname + ' written to disk.\n')

  # write PyMol input file
  fname = 'M' + str(monomers) + 'T' + str(timeval) + '.pml'
  fpymol = open(fname, 'w')
  fpymol.write('reinitialize\n')
  fpymol.write('viewport 640,480\n')
  for k in range(number_of_full_files + final_file):
    flnm = 'M' + str(monomers) + 'T' + str(timeval) + 'F' + str(k) + '.pdb'
    fpymol.write('load ' + flnm + '\n')
    fpymol.write('util.chainbow(\'' + flnm[0:len(flnm) - 4] + '\')\n')
    fpymol.write('show spheres, ' + flnm[0:len(flnm) - 4] + '\n')
  if(write_box == 'yes'):
    fpymol.write('load box.pdb\n')
    fpymol.write('color white, box\n')
    fpymol.write('show sticks, box\n')
  fpymol.write('reset\n')
  fpymol.close()
  sys.stdout.write('  Wrote out PyMol file ' + fname + '\n\n')

if(write_box == 'yes'):
  # write out box file
  scaled_L = sf * L
  fbox = open('box.pdb', 'w')
  fbox.write(pdb % ('HETATM', 1, 'C', 'ATM', 1, 0.0, 0.0, 0.0))
  fbox.write(pdb % ('HETATM', 2, 'C', 'ATM', 2, 0.0, 0.0, scaled_L))
  fbox.write(pdb % ('HETATM', 3, 'C', 'ATM', 3, scaled_L, scaled_L, 0.0))
  fbox.write(pdb % ('HETATM', 4, 'C', 'ATM', 4, scaled_L, scaled_L, scaled_L))
  fbox.write(pdb % ('HETATM', 5, 'C', 'ATM', 5, 0.0, scaled_L, 0.0))
  fbox.write(pdb % ('HETATM', 6, 'C', 'ATM', 6, 0.0, scaled_L, scaled_L)) 
  fbox.write(pdb % ('HETATM', 7, 'C', 'ATM', 7, scaled_L, 0.0, 0.0))
  fbox.write(pdb % ('HETATM', 8, 'C', 'ATM', 8, scaled_L, 0.0, scaled_L))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 1, 2))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 3, 4))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 1, 5))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 1, 7))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 2, 6))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 2, 8))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 3, 5))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 4, 6))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 4, 8))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 5, 6))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 7, 8))
  fbox.write('%6s%5d%5d\n' % ('CONECT', 3, 7))
  fbox.close()
  sys.stdout.write('box.pdb written to disk.\n\n')

if(len(blockfiles) == 1 and pymol_load == 'yes'):
  # view chains in PyMol
  cmd = 'env pymol ' + fname
  os.system(cmd)

"""This Python script computes the average number of original neighbors
   per chain at time t."""

import sys
import glob
import pypolymer

# time step (tau)
dt = 0.01

# generation of file names: (1) glob, (2) loop, (3) manual or (4) log scale
blockfiles_method = 1

if(blockfiles_method == 1):

  prefix = 'rings100a.'
  existingfiles = glob.glob(prefix + '*')
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]
  if 0:
    maxfiles = 30
    blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))]

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

  # number of files per decade
  n = 50

  # exponent of first and last decade (time steps)
  start_exp = 3
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
  sys.exit()

if(len(blockfiles) == 0):
  sys.stderr.write('ERROR: blockfiles is empty\n')
  sys.exit()

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                        n_N(t)                   \n')
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

# initialization
neighbors = []

xcm_zero_wrapped = []
ycm_zero_wrapped = []
zcm_zero_wrapped = []

xcm_zero = []
ycm_zero = []
zcm_zero = []

data = []

rgsq = {100:17.16, 200:30.77, 400:52.92, 800:87.61, 1600:145.61}
resq = {100:50.76, 200:88.75, 400:149.41, 800:242.41, 1600:401.73}

for file in blockfiles:

  xcm_t = []
  ycm_t = []
  zcm_t = []

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

  # additional parameters
  L = (total_particles / 0.85)**(1.0 / 3.0)
  M_rings27 = 27 * M_rings
  rcsq = resq[monomers]

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

  sys.stdout.write('Working on '+ file + '\n')

  # BEGIN: use first file to get initial values
  if(file == blockfiles[0]):

    # store timestep at t = 0
    time_zero = timestep
    if(timestep != 0):
      sys.stdout.write('WARNING: time is not 0 for initial file (' + str(timestep) +  ').\n\n')

    # compute COM of each chain at t = 0
    for i in range(M_rings):
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

      # store unwrapped COM positions at t = 0
      xcm_zero.append(xcm)
      ycm_zero.append(ycm)
      zcm_zero.append(zcm)

      # shift chain such that the central box is between [0, L]
      ix = int(xcm / L)
      iy = int(ycm / L)
      iz = int(zcm / L)
   
      # correct in case of negative positions
      if(xcm < 0.0): ix = ix - 1
      if(ycm < 0.0): iy = iy - 1
      if(zcm < 0.0): iz = iz - 1

      # wrap COM position
      xcm = xcm - ix * L
      ycm = ycm - iy * L
      zcm = zcm - iz * L

      # store central box COM positions at t = 0
      xcm_zero_wrapped.append(xcm)
      ycm_zero_wrapped.append(ycm)
      zcm_zero_wrapped.append(zcm)

    # translate central positions to 26 neighboring boxes
    for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          if(i**2 + j**2 + k**2 != 0):
	    for p in range(M_rings):
	      xcm_zero_wrapped.append(xcm_zero_wrapped[p] + i * L)
	      ycm_zero_wrapped.append(ycm_zero_wrapped[p] + j * L)
	      zcm_zero_wrapped.append(zcm_zero_wrapped[p] + k * L)
    if(len(xcm_zero_wrapped) != M_rings27): print 'Error: translation'

    # get neighbors at time zero
    for i in range(M_rings):
      tmp = []
      xcmi = xcm_zero_wrapped[i]
      ycmi = ycm_zero_wrapped[i]
      zcmi = zcm_zero_wrapped[i]
      for j in range(M_rings27):
        if(i != j):
          xcmj = xcm_zero_wrapped[j]
          ycmj = ycm_zero_wrapped[j]
          zcmj = zcm_zero_wrapped[j]
          rijsq = (xcmi - xcmj)**2 + (ycmi - ycmj)**2 + (zcmi - zcmj)**2
          if(rijsq <= rcsq):
            tmp.append(j)
      neighbors.append(tmp)
    if(len(neighbors) != M_rings): print 'Error: neighbors'
  # END: use first file to get initial values

  # compute COM of each chain
  for i in range(M_rings):
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

    # update positions of central box
    xcm_t.append(xcm_zero_wrapped[i] + xcm - xcm_zero[i])
    ycm_t.append(ycm_zero_wrapped[i] + ycm - ycm_zero[i])
    zcm_t.append(zcm_zero_wrapped[i] + zcm - zcm_zero[i])

  # translate central positions to 26 neighboring boxes
  for i in [-1, 0, 1]:
    for j in [-1, 0, 1]:
      for k in [-1, 0, 1]:
	if(i**2 + j**2 + k**2 != 0):
	  for p in range(M_rings):
	    xcm_t.append(xcm_t[p] + i * L)
	    ycm_t.append(ycm_t[p] + j * L)
	    zcm_t.append(zcm_t[p] + k * L)

  # count neighbors at time t
  ct = 0
  for i in range(M_rings):
    xcmi = xcm_t[i]
    ycmi = ycm_t[i]
    zcmi = zcm_t[i]
    for j in neighbors[i]:
	xcmj = xcm_t[j]
	ycmj = ycm_t[j]
	zcmj = zcm_t[j]
	rijsq = (xcmi - xcmj)**2 + (ycmi - ycmj)**2 + (zcmi - zcmj)**2
	if(rijsq <= rcsq):
          ct = ct + 1
  ave_original_neighbors_per_chain = ct / float(M_rings)
  data.append((dt * (timestep - time_zero), ave_original_neighbors_per_chain))

f = open('orig_neighbors' + str(monomers) + '.dat', 'w')
f.write('# time (tau)  ave_original_neighbors_per_chain\n')
for rec in data:
  f.write('%.2f %.3f\n' % (rec[0], rec[1])) 
f.close()

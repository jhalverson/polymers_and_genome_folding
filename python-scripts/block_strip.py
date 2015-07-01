"""This script prepares block files for the Fortran
   code that computes the self-excluded RDF."""

chains = 200

istr = 62600000
iend = 66000000
incr = 100000

times = []; blockfiles = []
for i in range(istr, iend + incr, incr):
  blockfiles.append('RG_N200M200a.0' + str(i))
  times.append(str(i))

# loop over block files
for file, time in zip(blockfiles, times):

  f = open(file)

  # find line where total particles is found
  line = ''
  while not '{n_part ' in line:
    line = f.readline()

  sline = line.replace('}','')
  total_particles = int(sline.split()[1])
  monomers = total_particles / chains
  
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

  flnm = str(time) + '.pos'
  F = open(flnm, 'w')
  for i in range(total_particles):
    outline = '%g %g %g \n' % (x[i], y[i], z[i])
    F.write(outline)

  F.close()

  print 'monomers=', monomers, '  particles=', total_particles, '  filename=', flnm

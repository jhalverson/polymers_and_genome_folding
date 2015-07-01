import sys
import os
import glob

files = glob.glob('dump_conf*')

for file in files:

  os.system('wc -l ' + file + ' > lines.dat')
  f = open('lines.dat')
  lines = int(f.readline().split()[0])
  f.close()

  f = open(file)
  line = ''
  while not 'ITEM: NUMBER OF ATOMS' in line:
    line = f.readline()
  total_particles = int(f.readline())
  f.close()

  f = open(file)
  configs = int(lines / (total_particles + 9))
  ct = 0
  while (ct < configs):

    line = ''
    while not 'ITEM: TIMESTEP' in line:
      line = f.readline()
    timestep = int(f.readline())

    line = ''
    while not 'ITEM: NUMBER OF ATOMS' in line:
      line = f.readline()
    total_particles = int(f.readline())

    line = ''
    while not 'ITEM: BOX BOUNDS' in line:
      line = f.readline()
    xmin, xmax = map(float, f.readline().split())
    L = xmax - xmin

    if(abs(1.0 - (total_particles / L**3) / 0.85) > 0.01):
      print 'density is different: ', total_particles / L**3
      sys.exit()

    line = ''
    while not 'ITEM: ATOMS' in line:
      line = f.readline()
    coords = []
    for i in range(total_particles):
      sline = f.readline()
      p = int(sline.split()[0])
      x = float(sline.split()[2]) + L * int(sline.split()[5])
      y = float(sline.split()[3]) + L * int(sline.split()[6])
      z = float(sline.split()[4]) + L * int(sline.split()[7])
      coords.append((p, x, y, z))

    coords.sort(lambda u, v: cmp(u[0], v[0]))

    print ct,'of', configs, timestep, total_particles, L, coords[0], coords[-1]
    g = open(str(timestep) + '.pos', 'w')
    for i in range(total_particles):
      g.write('%10.5f %10.5f %10.5f\n' % (coords[i][1], coords[i][2], coords[i][3]))
    g.close()

    ct = ct + 1

  f.close()

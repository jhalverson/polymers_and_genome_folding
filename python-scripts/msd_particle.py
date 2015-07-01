import sys
import os
import glob

# ignore data before this time step
ignore = 1000000

# time step (tau)
delta_t = 0.01

files = glob.glob('*RADIUS*')
files = filter(lambda u: u[u.rindex('.') + 1:].isdigit() and \
                         'restart' not in u and \
                         int(u[u.rindex('.') + 1:]) >= ignore, files)
files.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))

txyz = []
for file in files:
  f = open(file)
  data = f.readlines()
  f.close()

  timestep = data[1].strip()
  particle_id = data[3].strip()
  for d in data:
    s = d.split()
    if particle_id == s[0] and len(s) == 4:
      txyz.append((int(timestep), float(s[1]), float(s[2]), float(s[3])))

ct = {}
g3 = {}
for i, (t0, x0, y0, z0) in enumerate(txyz):
  for t, x, y, z in txyz[i:]:
    dt = t - t0
    if not ct.has_key(dt):
      ct[dt] = 0
      g3[dt] = 0.0
    ct[dt] += 1
    g3[dt] += (x - x0)**2 + (y - y0)**2 + (z - z0)**2

f = open('msd_' + files[0][:files[0].rindex('.')]+ '.dat', 'w')
f.write('# time (tau)  g3 / sigma**2\n')
for dt in sorted(ct.keys()):
  f.write('%8d %9.1f\n' % (delta_t * dt, g3[dt] / ct[dt]))
f.close()

"""This python script combines stress data files from each job step
   into a single file per component."""

import sys
import glob

# which contributions
contrib = ['ke', 'pair', 'bond', 'angle', 'total']

# which components
comp = ['pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz']

# job name
jobname = 'rings400a'

# shift time (time steps)
shift = 500

# dictionary of stress tensor components (do not change this)
press = {'pxx': 0, 'pyy': 1, 'pzz': 2, 'pxy': 3, 'pxz': 4, 'pyz': 5}

tlist = []
plist = []
for c in contrib:
  files = glob.glob('*' + c + '*')
  files.sort(lambda u, v: cmp(int(u[u.index('dat.') + 4:]), \
                              int(v[v.index('dat.') + 4:])))
  for fl in files:
    f = open(fl)
    lines = f.readlines()
    f.close()
    lines = filter(lambda u: not u.strip().startswith('#'), lines)
    for i in range(0, len(lines), 7):
      timestep, num = lines[i].split()
      if(c == contrib[0]): tlist.append(timestep)
      num, pxx = lines[i + 1].split()
      num, pyy = lines[i + 2].split()
      num, pzz = lines[i + 3].split()
      num, pxy = lines[i + 4].split()
      num, pxz = lines[i + 5].split()
      num, pyz = lines[i + 6].split()
      pij = [pxx, pyy, pzz, pxy, pxz, pyz]
      for d in comp:
        plist.append(pij[press[d]])

# have the timestep list and one long list for the stress components

if(len(plist) != len(tlist) * len(comp) * len(contrib)):
  sys.stderr.write('ERROR: Length of list is wrong.\n')

for k, d in enumerate(comp):
  outfile = jobname + '_stress.' + d
  f = open(outfile, 'w')
  f.write('# timestep ' + ' '.join(contrib) + '\n')
  for i in range(len(tlist)):
    timestep = str(int(tlist[i]) - shift)
    s = ''
    for j in range(len(contrib)):
      p = j * len(tlist) * len(comp) + i * len(comp) + k
      s = s + plist[p] + ' '
    f.write(timestep + ' ' + s + '\n')
  f.close()
  sys.stdout.write(outfile + ' has been written to disk.\n')

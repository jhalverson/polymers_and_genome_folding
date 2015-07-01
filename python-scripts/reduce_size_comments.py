"""This Python script extracts the off-diagonal components
   from the data of Grest."""

import sys

first_reset_step = 14
laststep = -1e99

pfiles = []
for i in range(1, 12):
  pfiles.append('stress.dat.' + str(i))
print pfiles

p = []
for pf in pfiles:

  f = open(pf)
  data_c = f.readlines()
  f.close()

  data = filter(lambda u: not u.strip().startswith('#'), data_c)

  step = int(pf[pf.find('t.') + 2:])
  print step

  for i in range(len(data)):
    if(step >= first_reset_step):
      timestep = int(data[i].split()[0]) + laststep
      s = str(timestep) + ' ' + data[i].split()[4] + ' ' + data[i].split()[5] + ' ' + data[i].split()[6] + '\n'
    else:
      s = data[i].split()[0] + ' ' + data[i].split()[4] + ' ' + data[i].split()[5] + ' ' + data[i].split()[6] + '\n'
    p.append(s)

f = open('200_1_pxy_pxz_pyz.dat', 'w')
for i in range(len(p)):
  f.write(p[i])
f.close()

for i in range(len(p)):
  if(100000000 + (i+1)*1000 != int(p[i].split()[0])):
    print i, p[i]
    sys.exit()

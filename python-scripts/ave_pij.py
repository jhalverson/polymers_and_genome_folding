"""This Python script reads in the large stress data files and
   reduces the data by applying a time-averaging scheme."""

import glob

# averaging window (MD steps)
ave_window = 1000

# job name
jobname = '400a_every'

pfiles = glob.glob('../../pij.dat.*')
pfiles.sort(lambda u, v: cmp(int(u[u.index('dat.') + 4:]), int(v[v.index('dat.') + 4:])))
print pfiles

iwindow = 0
for pf in pfiles:

  f = open(pf)
  data_c = f.readlines()
  f.close()

  jobstep = int(pf[pf.index('dat.') + 4:])

  print pf
  print data_c[3], (jobstep - 2) * 7e7 + 5e7 + 1
  print data_c[-1], (jobstep - 1) * 7e7 + 5e7

  sz = len(data_c)
  number_of_windows = (sz - 3) / ave_window
  print (sz - 3), number_of_windows
  print '\n'

  f = open('ave' + str(ave_window) + '_' + jobname + '.dat.' + str(jobstep), 'a')

  for i in xrange(number_of_windows):
    ave_pxy = 0.0
    ave_pxz = 0.0
    ave_pyz = 0.0
    for j in xrange(ave_window):
      p = i * ave_window + j + 3
      ave_pxy = ave_pxy + float(data_c[p].split()[1])
      ave_pxz = ave_pxz + float(data_c[p].split()[2])
      ave_pyz = ave_pyz + float(data_c[p].split()[3])
    ave_pxy = ave_pxy / ave_window
    ave_pxz = ave_pxz / ave_window
    ave_pyz = ave_pyz / ave_window
    iwindow = iwindow + 1
    f.write('%d %.7e %.7e %.7e\n' % (iwindow * ave_window, ave_pxy, ave_pxz, ave_pyz))
  f.close()

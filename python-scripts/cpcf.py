"""This Python script computes a cross product correlation
   function averaged over all the chains in the system."""

import sys

# data file
fin = 'ab100c_n5.dat'

# spacing (tau)
spacing = 5000

# number of chains
chains = 200

# number of files
files = 4000

# maximum correlation time
t_cor = 2000

f = open(fin)
line_first = f.readline()
remove_comment = False
if(line_first.strip().startswith('#')):
  remove_comment = True
f.close()

c = []
f = open(fin)
if(remove_comment): f.readline()
for i in range(files):
  c_per_file = []
  for j in range(chains):
    line = f.readline().split()
    a = (float(line[2]), float(line[3]), float(line[4]))
    b = (float(line[5]), float(line[6]), float(line[7]))
    v = (float(line[8]), float(line[9]), float(line[10]))
    c_per_file.append(v)
  c.append(c_per_file)
f.close()

# check
if(len(c) != files): sys.stdout.write('Length of c is wrong.\n')

# total number of data values
t_run = files

# compute time-correlation function of each chain using the direct method
cf_all = []
for i in range(chains):
  cf = []
  for t in range(t_cor):
    t_max = t_run - t
    sum_cc = 0.0
    for t0 in range(t_max):
      sum_cc = sum_cc + c[t0][i][0] * c[t0 + t][i][0] \
                      + c[t0][i][1] * c[t0 + t][i][1] \
                      + c[t0][i][2] * c[t0 + t][i][2]
    cf.append(sum_cc / t_max)
  cf_all.append(cf)

# average results
cf_ave = []
for t in range(t_cor):
  cc_ave = 0.0
  for i in range(chains):
    cc_ave = cc_ave + cf_all[i][t]
  cf_ave.append(cc_ave / chains)

# write data to file
outfile = fin + '.acf'
f = open(outfile, 'w')
for t in range(t_cor):
  f.write('%8d %.6f\n' % (spacing * t, cf_ave[t] / cf_ave[0]))
f.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

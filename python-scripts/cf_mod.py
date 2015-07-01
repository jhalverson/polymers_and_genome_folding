"""This script computes the time correlation function."""

# use sys for printing compatibility between Python 2 & 3
import sys

# data file
fin = 'L23_all_800a.dat'

# chains
chains = 200

# read and store file contents
f = open(fin)
lines = f.readlines()
f.close()

# store number of lines in file
lines_len_init = len(lines)

# remove comment lines
lines = filter(lambda u: not u.strip().startswith('#'), lines)

# display number of lines and number removed
sys.stdout.write('Number of lines: ' + str(len(lines)) + '\n')
num_comments = lines_len_init - len(lines)
if(num_comments != 0): sys.stdout.write('Number of comment lines removed: ' + str(num_comments) + '\n')

# load data
a = []
for j in range(len(lines)):
  line = lines[j]
  s = line.split()
  if(j == 0): slen = len(s)
  if(len(s) != slen): sys.stderr.write('WARNING: number of columns has changed.\n')
  # ignore time columns
  s = map(lambda u: float(u), s[2:])
  if(len(s) != chains): sys.stderr.write('WARNING: number of columns not equal to chains.\n')
  a.append(s)

# total number of data values per chain
t_run = len(lines)

# maximum correlation time (in units of configuration files)
t_cor = t_run / 2
sys.stdout.write('t_run, t_cor, t_run / t_cor = %d %d %d\n' % (t_run, t_cor, t_run / t_cor))

# time difference between files (in units of tau)
t_diff = 5000

# get average values and variance per chain
a_ave = []
a_var = []
for i in range(chains):
  # average
  a_sum = 0.0
  for row in a:
    a_sum = a_sum + row[i]
  a_ave.append(a_sum / len(a))
  # variance
  a_diff_sum = 0.0
  for row in a:
    a_diff_sum = a_diff_sum + (row[i] - a_ave[i])**2
  a_var.append(a_diff_sum / len(a))

# compute time correlation function of each chain using the direct method
cf_all = []
for i in range(chains):
  cf = []
  for t in range(t_cor):
    t_max = t_run - t
    sum_aa = 0.0
    for t0 in range(t_max):
      sum_aa = sum_aa + (a[t0][i] - a_ave[i]) * (a[t0 + t][i] - a_ave[i])
    cf.append(sum_aa / t_max)
  cf_all.append(cf)

# average single chain acf's into a master acf
cf_ave = []
for t in range(t_cor):
  aa_ave = 0.0
  for i in range(chains):
    aa_ave = aa_ave + cf_all[i][t]
  cf_ave.append(aa_ave / chains)

# check that normalization is close to correct value
sum_var = 0.0
for i in range(chains):
  sum_var = sum_var + a_var[i]
sum_var = sum_var / chains
if(abs(sum_var / cf_ave[0] - 1.0) > 0.01):
  sys.stderr.write('WARNING: Problem with normalization.\n')

# write data to file
outfile = 'acf_' + fin[:]
f = open(outfile, 'w')
for t in range(t_cor):
  f.write('%10d %10.4f\n' % (t_diff * t, cf_ave[t] / cf_ave[0]))
f.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

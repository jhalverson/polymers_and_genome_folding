"""This script computes the time correlation function
   for the indicated columnar data."""

# use sys for printing compatibility between Python 2 & 3
import sys

# data file
fin = 'eigenvalues_cf1600a.dat'

# column 1
c1 = 0

# column 2
c2 = 4

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
for i in range(len(lines)):
  line = lines[i]
  s = line.split()
  if(i == 0): slen = len(s)
  if(len(s) != slen): sys.stderr.write('WARNING: number of columns has changed.\n')
  a.append(float(s[c2]))

# normalization
Q = sum(map(lambda u: u**2, a)) / len(a)

# total number of data values
t_run = len(lines)

# maximum correlation time
t_cor = t_run / 1
sys.stdout.write('total length divided by expected correlation length: %.1f\n' % (float(t_run) / t_cor))

# compute time-correlation function using the direct method
cf = []
for t in range(t_cor):
  t_max = t_run - t
  sum_aa = 0.0
  for t0 in range(t_max):
    sum_aa = sum_aa + a[t0] * a[t0 + t]
  cf.append((t, sum_aa / t_max / Q))

# write data to file
outfile = 'acf_L13_1.dat'
f = open(outfile, 'w')
for i in range(len(cf)):
  f.write('%10d %10.4f\n' % (cf[i][0], cf[i][1]))
f.close()
sys.stdout.write(outfile + ' has been written to disk.\n')

"""This script reads in a g123 file and outputs n points
   for each decade on a logarithmic scale. It may be necessary
   to remove the first and last point of the new data set."""

import sys
import os
import math

# number of points to keep per decade
n = 10

# set the base
b = math.e

# load data into Gnuplot ('yes' or 'no')
gnuplot_load = 'yes'

# generate file name: (1) glob or (2) manual
fname_method = 1

if(fname_method == 1):

  import glob
  g123file = glob.glob('*.g123')

  # check that only one file was found
  if(len(g123file) != 1):
    sys.stdout.write('ERROR: zero or more than one g123 file was found.\n\n')
    sys.exit(1)

  # convert list to string
  g123file = g123file[0]
  sys.stdout.write('\n\nOperating on ' + g123file + ' ...\n')

elif(fname_method == 2):

  # manually enter name of g123 file on next line
  g123file = 'RG_N200M1600.g123'
  sys.stdout.write('\n\nOperating on ' + g123file + ' ...\n')
 
else:

  sys.stdout.write('ERROR: fname_method\n')
  sys.exit(1)

# store data in g123
g123 = []
f = open(g123file, 'r')
for line in f:
  g123.append(line.split())
f.close()

# show the user the first and last line of g123file
sys.stdout.write('  Each column of ' + g123file + ' has ' + str(len(g123)) + ' values.\n')
sys.stdout.write('  First line: ' + ' '.join(g123[0]) + '\n')
sys.stdout.write('  Last line:  ' + ' '.join(g123[len(g123) - 1]) + '\n')

# starting decade
start10 = '%e' % float(g123[0][0])
start_exp = int(start10[start10.index('e') + 1:]) - 1

# ending decade
end10 = '%e' % float(g123[-1][0])
end_exp = int(end10[end10.index('e') + 1:]) + 1

# create list of numbers uniformly spaced on a log scale from 10^start_exp to 10^end_exp
uniform_time = []
number_of_decades = end_exp - start_exp
start = math.log(10**start_exp)
end = math.log(10**end_exp)
step = (end - start) / (n * number_of_decades - 1)
for i in range(n * number_of_decades):
  uniform_time.append(b**(start + step * i))

# loop over uniformly-spaced values and find closest value in g123
uniform_g123_raw = []
for t in uniform_time:
  mindiff = float('inf')
  for line in g123:
    diff = math.fabs(t - float(line[0]))
    if(diff < mindiff):
      closest_line = line
      mindiff = diff
  uniform_g123_raw.append(closest_line)

# remove duplicates
uniform_g123 = [uniform_g123_raw[0]]
for line in uniform_g123_raw[1:]:
  if(line != uniform_g123[-1]):
    uniform_g123.append(line)

# write output file
outfile = g123file + '.out'
fout = open(outfile, 'w')
for line in uniform_g123:
  fout.write(' '.join(line) + '\n')
fout.close()
sys.stdout.write('  Output file ' + outfile + ' has been written to disk.\n')

# write a Gnuplot input file
pltfile = g123file + '.plt'
f = open(pltfile, 'w')
f.write('set xlabel "TIME (TAU)"\n')
f.write('set ylabel "MSD (SIGMA * SIGMA)"\n')
f.write('set key left\n')
f.write('set grid\n')
f.write('set log xy\n')
f.write('set xrange[1e3: 1e7]\n')
f.write('set yrange[0.1: 1000]\n')
f.write('set pointsize 2\n')
f.write('set mxtics 10\n')
f.write('set mytics 10\n')
f.write('plot \\\n')
f.write('"' + g123file + '" using 1:2 with points title "g1", \\\n')
f.write('"' + g123file + '" using 1:3 with points title "g2", \\\n')
f.write('"' + g123file + '" using 1:4 with points title "g3", \\\n')
f.write('"' + outfile  + '" using 1:2 with points title "g1", \\\n')
f.write('"' + outfile  + '" using 1:3 with points title "g2", \\\n')
f.write('"' + outfile  + '" using 1:4 with points title "g3", \\\n')
f.write('(x/10000) title "t", (x/10)**0.5 title "t^0.5"')
f.close()
sys.stdout.write('  Gnuplot input file ' + pltfile + ' has been written to disk.\n\n')

if(gnuplot_load == 'yes'):
  cmd = 'gnuplot -persist ' + pltfile
  os.system(cmd)

"""This Python script compute the average and standard deviation
   of data written to a LAMMPS log file."""

import sys
import numpy as np

f = open('log.lammps')

line = ''
while not 'Step Temp Press PotEng' in line:
  line = f.readline()
sys.stdout.write('\n\n' + line + '\n')

# initialize arrays
T = []; P = []; E = []

line = f.readline()
while line.split()[0].isdigit():
  stp = int(line.split()[0])
  T.append(float(line.split()[1]))
  P.append(float(line.split()[2]))
  E.append(float(line.split()[3]))
  sys.stdout.write('(step, T, P, E) = (%d, %.3f, %.3f, %.3f)\n' % (stp, T[-1], P[-1], E[-1]))
  line = f.readline()
  if(line == ''): break
f.close()

T = np.array(T)
P = np.array(P)
E = np.array(E)

ave_T = np.mean(T); std_T = np.std(T)
ave_P = np.mean(P); std_P = np.std(P)
ave_E = np.mean(E); std_E = np.std(E)

sys.stdout.write('\n(T, P, E) = %.2f (%.3f)  %.2f (%.3f)  %.2f (%.3f)\n' % \
                (ave_T, std_T, ave_P, std_P, ave_E, std_E))

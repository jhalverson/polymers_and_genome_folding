import os
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pypolymer

plt.rcParams['text.usetex'] = True
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['font.size'] = 8
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['legend.fontsize'] = 6
plt.rcParams['legend.borderaxespad'] = 1
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['lines.markersize'] = 3

w = 8.0 / 2.54
h = 0.2

for i in range(0, 360000 + 1000, 1000):
#for i in range(0, 1000, 1000):
  fig = plt.figure(i + 1, figsize=(w, h))
  plt.figtext(0.05, 0.15, r'$N=200,~M_{\mathrm{rings}}=200,~M_{\mathrm{linear}}=26$')
  plt.figtext(0.78, 0.33, r'$t=' + str(i/100) + r'~\tau~~~~~~~~~~~1$')

  outfile = 'labels' + str(i)
  plt.savefig(outfile + '.png', dpi=203.5)
#  os.system('xv ' + outfile + '.png&\n')

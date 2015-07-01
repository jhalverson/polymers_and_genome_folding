import os
import glob
import time
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import sys
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

def half(a_):
  h = []
  if(len(a_) % 2 == 0):
    for i in range(0, len(a_), 2):
      h.append(0.5 * (a_[i] + a_[i+1]))
  else:
    for i in range(0, len(a_) - 1, 2):
      h.append(0.5 * (a_[i] + a_[i + 1]))
  return h

def block(b_, m_):
  for i in range(0, m_):
    b_ = half(b_)
  return b_

if 1:

  gamma_dot = 1e-6
  files = glob.glob('stress.dat.*')
  files.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  print files
  pxz = []
  for file in files:
    if('norm' in file):
      pxz_ = pypolymer.get_data(file, 5)
    else:
      pxz_ = pypolymer.get_data(file, 2)
    pxz.extend(pxz_)
  eta = np.array(pxz) / -gamma_dot
  print 'ave=', np.mean(eta), 'std=', np.std(eta)

  f = open('blocks.dat', 'w')
  ct = 0
  sg = []
  sg_err = []
  for M in range(50):
    z = block(eta, M)
    Lp = len(z)
    if(Lp < 2): break
    ave = np.mean(z)
    var = np.var(z)
    std = (var / (Lp - 1))**0.5
    std2 = std * (1.0 / (2.0 * (Lp - 1))**0.5)
    sg.append(std)
    sg_err.append(std2)
    print M, Lp, ave, std, std2
    f.write('%d %d %g %g %g\n' % (M, Lp, ave, std, std2))
    ct = ct + 1
  f.close()

w = 8.0 / 2.54
h = (3.0 / 4.0) * w

fig = plt.figure(1, figsize=(w, h))
ax = plt.subplot(1, 1, 1)
plt.errorbar(range(ct), sg, yerr=sg_err, fmt='ko')
plt.title(r'$\mathrm{Error~analysis~for~Rings}~N=800,~\dot\gamma=10^{-6}$')
plt.xlabel(r'$M$')
plt.ylabel(r'$\sigma_{\eta}$')
plt.figtext(0.02, 0.95, os.getcwd().replace('_', '\_'), fontsize=4, color='0.75')
plt.figtext(0.02, 0.93, time.asctime(), fontsize=4, color='0.75')

fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.8, wspace=0.0, hspace=0.0)

for ax in fig.axes:
  plt.setp(ax.patch, color=(0.9, 0.9, 0.9))

plt.savefig('blocks.png', dpi=200)
os.system('xv blocks.png&\n')

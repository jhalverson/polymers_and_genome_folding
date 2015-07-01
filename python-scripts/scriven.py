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

keep_last_percent = 50 # percentage to keep starting from the end
gamma_dot = 1e-4

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

# ignore early data
start = len(pxz) * (1.0 - float(keep_last_percent) / 100)
pxz = pxz[int(start):]

eta = np.array(pxz) / -gamma_dot
print 'ave=', np.mean(eta), 'std=', np.std(eta)

pxz = np.array(pxz) / -gamma_dot
pxz = pxz.tolist()

block_pre = filter(lambda u: u <= len(pxz), [2**i for i in range(25)])
if(block_pre[-1] != len(pxz)): block_pre.append(len(pxz))

ave_plot = []
block_plot = []

f = open('blocks.dat', 'w')
for block in block_pre:
  blocks = len(pxz) / block
  pxz_blocked = []
  for b in range(blocks):
    ave = 0.0
    for i in range(block):
      ave = ave + pxz[b * block + i]
    pxz_blocked.append(ave / block)
  ave_plot.append(np.std(np.array(pxz_blocked)))
  block_plot.append(block)
  f.write('%d %d %.4e\n' % (block, blocks, np.std(np.array(pxz_blocked))))
f.close()

w = 8.0 / 2.54
h = (3.0 / 4.0) * w

fig = plt.figure(1, figsize=(w, h))
ax = plt.subplot(1, 1, 1)
plt.semilogx(block_plot, ave_plot, 'ko', basex=2)
plt.title(r'$\mathrm{Error~analysis~for~Rings}~N=100,~\dot\gamma=10^{-4}$')
plt.xlabel(r'$\mathrm{block~size}$')
plt.ylabel(r'$\sigma_{\eta}$')
plt.figtext(0.02, 0.95, os.getcwd().replace('_', '\_'), fontsize=4, color='0.75')
plt.figtext(0.02, 0.93, time.asctime(), fontsize=4, color='0.75')

fig.subplots_adjust(left=0.2, right=0.85, bottom=0.2, top=0.8, wspace=0.0, hspace=0.0)

for ax in fig.axes:
  plt.setp(ax.patch, color=(0.9, 0.9, 0.9))

plt.savefig('blocks.png', dpi=200)
os.system('xv blocks.png&\n')

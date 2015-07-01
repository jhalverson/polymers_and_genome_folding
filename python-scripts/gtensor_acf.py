"""This Python script computes the normalized autocorrelation function
   of Rg2 and the ratios of the eigenvalues of the gyration tensor."""

import numpy as np

def autocorr(x):
  """Takes a list and returns the autocorrelation function as a numpy array."""
  s = np.fft.fft(x, n = 2 * len(x))
  acf = np.real(np.fft.ifft(s * np.conjugate(s)))[0:len(x)]
  for i in range(acf.size):
    acf[i] = acf[i] / (acf.size - i)
  return acf

qty = 'L23'
#qty = 'L13'
#qty = 'RGsq'
monomers = 1600
chains = 200
single_ave = False

f = open(qty + '_all_' + str(monomers) + 'a.dat')
data = f.readlines()
f.close()
numlines = len(data)
interval = int(data[1].split()[0])
print 'numlines =', numlines
print 'interval =', interval
print 'monomers =', monomers

rg2 = []
for i in range(chains):
  tmp = []
  for j in range(numlines):
    tmp.append(float(data[j].split()[i + 2]))
  rg2.append(tmp)

del data

ave_rg2 = 0.0
for i in range(chains):
  ave_rg2 += np.mean(np.array(rg2[i]))
ave_rg2 = ave_rg2 / chains
print 'ave =', ave_rg2

for i in range(chains):
  tmp = np.mean(np.array(rg2[i]))
  for j in range(numlines):
    if(single_ave):
      rg2[i][j] = rg2[i][j] - tmp
    else:
      rg2[i][j] = rg2[i][j] - ave_rg2

acf_rg2 = np.zeros(numlines, dtype = np.float64)
for i in range(chains):
  acf_rg2 += autocorr(rg2[i])
acf_rg2 = acf_rg2 / chains
acf_rg2 = acf_rg2 / acf_rg2[0]

f = open(qty + '_acf' + str(monomers) + '.dat', 'w')
if(single_ave): f = open(qty + '_acf' + str(monomers) + '_single.dat', 'w')
f.write('# time (tau)  acf\n')
for i in range(acf_rg2.size):
  f.write('%d %8.4f\n' % (i * interval, acf_rg2[i]))
f.close()

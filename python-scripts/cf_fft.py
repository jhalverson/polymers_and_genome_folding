"""This Python script computes the time correlation function
   using the fast Fourier transform method. See Allen & Tildesly
   for an explanation. There are times to use the direct approach
   and times to use the FFT method."""

# use sys for printing compatibility between Python 2 & 3
import sys
import numpy as np
import os
import time

# data file
fin = 'data.dat'

# read and store file contents
f = open(fin)
lines = f.readlines()
f.close()

# load data
a = []
for i in range(len(lines)):
  line = lines[i]
  s = line.split()
  a.append(float(s[1]))

def autocorr(x):
  """Takes a list and returns the autocorrelation function as a numpy array."""
  s = np.fft.fft(x, n = 2 * len(x))
  acf = np.real(np.fft.ifft(s * np.conjugate(s)))[0:len(x)]
  for i in range(acf.size):
    acf[i] = acf[i] / (acf.size - i)
  return acf

b = autocorr(a)
print b

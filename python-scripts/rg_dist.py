"""This Python script computes P(Rg), P(Re), P(Rg2) and P(Re2)
   for linear chains or rings."""

# use sys for printing compatibility between Python 2 & 3
import sys
import math
import numpy as np
import pypolymer

# number of chains M
chains = 200

# rings or linear
topology = 'rings'

# run
case = ''

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  import glob
  blockfiles = glob.glob('rings400a.*')
  blockfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), blockfiles)
  blockfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))

  # work with a subset
  blockfiles = blockfiles[:]

elif(blockfiles_method == 2):

  istr = 1000000
  iend = 2000000000
  incr = 1000000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('../sqt/stripped/' + str(i) + '.pos')

elif(blockfiles_method == 3):

  blockfiles = ['rings400a.5000000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit()

sys.stdout.write('Total number of files: ' + str(len(blockfiles)) + '\n\n')

# initialize histogram for RG
dRG = 0.5
upper_bound_RG = 30.0
ibinmax = int(upper_bound_RG / dRG)
histRG = []
for i in range(ibinmax):
  histRG.append(0)

# initialize histogram RE
dRE = 1.0
upper_bound_RE = 30.0
ibinmax = int(upper_bound_RE / dRE)
histRE = []
for i in range(ibinmax):
  histRE.append(0)

# initialize histogram for RG2
dRG2 = 5.0
upper_bound_RG2 = 500.0
ibinmax = int(upper_bound_RG2 / dRG2)
histRG2 = []
for i in range(ibinmax):
  histRG2.append(0)

# initialize histogram for RE2
dRE2 = 10.0
upper_bound_RE2 = 500.0
ibinmax = int(upper_bound_RE2 / dRE2)
histRE2 = []
for i in range(ibinmax):
  histRE2.append(0)

# loop over block files
for file in blockfiles:

  fmanager = pypolymer.read_blockfile(file)
  timestep = fmanager[0]
  total_particles = fmanager[1]
  M_rings = fmanager[2]; M_linear = fmanager[3]
  monomers = fmanager[4]; monomers_half = fmanager[5]
  x = fmanager[6]; y = fmanager[7]; z = fmanager[8]

  if(file == blockfiles[0]):
    total_particles1st = total_particles
    sys.stdout.write('System parameters:\n')
    sys.stdout.write('  M_rings = ' + str(M_rings) + '\n')
    sys.stdout.write('  M_linear = ' + str(M_linear) + '\n')
    sys.stdout.write('  N = ' + str(monomers) + '\n')
    sys.stdout.write('  N/2 = ' + str(monomers_half) + '\n\n')
  if(total_particles1st != total_particles):
    sys.stderr.write('Change in total_particles (' + str(total_particles) + ').\n')
    sys.exit()

  sys.stdout.write('Working on ' + file + '\n')

  # begin calculation
  for i in range(chains):

    # compute center-of-mass
    xcm = 0.0; ycm = 0.0; zcm = 0.0
    for j in range(monomers):
      p = i * monomers + j
      xcm = xcm + x[p]
      ycm = ycm + y[p]
      zcm = zcm + z[p]
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers

    # RG and RG2
    radsq = 0.0
    for j in range(monomers):
      p = i * monomers + j
      radsq = radsq + (x[p] - xcm)**2 + (y[p] - ycm)**2 + (z[p] - zcm)**2
    RG2 = radsq / monomers
    RG = math.sqrt(RG2)

    # RE and RE2
    spansq = 0.0
    if(topology == 'rings'):
      # compute spanning distances squared
      for j in range(monomers_half):
	p1 = i * monomers + j
	p2 = p1 + monomers_half
	spansq = spansq + (x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1] - z[p2])**2
    else:
      # compute end-to-end distances squared
      p1 = i * monomers
      p2 = p1 + monomers - 1
      spansq = spansq + (x[p1] - x[p2])**2 + (y[p1] - y[p2])**2 + (z[p1] - z[p2])**2
    if(topology == 'rings'):
      # compute diameter squared
      RE2 = spansq / monomers_half
      RE = math.sqrt(RE2)
    else:
      # compute end-to-end distance squared
      RE2 = spansq
      RE = math.sqrt(RE2)

    # RG
    ibin = int(RG / dRG)
    if(ibin < len(histRG)):
      histRG[ibin] = histRG[ibin] + 1
    else:
      sys.stderr.write('ERROR: ibinmax too small ' + str(RG))

    # RE
    ibin = int(RE / dRE)
    if(ibin < len(histRE)):
      histRE[ibin] = histRE[ibin] + 1
    else:
      sys.stderr.write('ERROR: ibinmax too small ' + str(RE))

    # RG2
    ibin = int(RG2 / dRG2)
    if(ibin < len(histRG2)):
      histRG2[ibin] = histRG2[ibin] + 1
    else:
      sys.stderr.write('ERROR: ibinmax too small ' + str(RG2))

    # RE2
    ibin = int(RE2 / dRE2)
    if(ibin < len(histRE2)):
      histRE2[ibin] = histRE2[ibin] + 1
    else:
      sys.stderr.write('ERROR: ibinmax too small ' + str(RE2))

# RG
pdfRG = np.array(histRG) / float(chains * len(blockfiles))
area = np.sum(pdfRG * dRG)
pdfRG = pdfRG / area

outfile = 'rg_dist' + str(monomers) + case + '.dat'
fout = open(outfile, 'w')
for i in range(pdfRG.size):
  minRG = i * dRG
  maxRG = minRG + dRG
  aveRG = minRG + dRG / 2.0
  outline = '%g %g %g %g\n' % (minRG, maxRG, aveRG, pdfRG[i])
  fout.write(outline)
fout.close()

# RE
pdfRE = np.array(histRE) / float(chains * len(blockfiles))
area = np.sum(pdfRE * dRE)
pdfRE = pdfRE / area

outfile = 're_dist' + str(monomers) + case + '.dat'
fout = open(outfile, 'w')
for i in range(pdfRE.size):
  minRE = i * dRE
  maxRE = minRE + dRE
  aveRE = minRE + dRE / 2.0
  outline = '%g %g %g %g\n' % (minRE, maxRE, aveRE, pdfRE[i])
  fout.write(outline)
fout.close()

# RG2
pdfRG2 = np.array(histRG2) / float(chains * len(blockfiles))
area = np.sum(pdfRG2 * dRG2)
pdfRG2 = pdfRG2 / area

outfile = 'rg2_dist' + str(monomers) + case + '.dat'
fout = open(outfile, 'w')
for i in range(pdfRG2):
  minRG2 = i * dRG2
  maxRG2 = minRG2 + dRG2
  aveRG2 = minRG2 + dRG2 / 2.0
  outline = '%g %g %g %g\n' % (minRG2, maxRG2, aveRG2, pdfRG2[i])
  fout.write(outline)
fout.close()

# RE2
pdfRE2 = np.array(histRE2) / float(chains * len(blockfiles))
area = np.sum(pdfRE2 * dRE2)
pdfRE2 = pdfRE2 / area

outfile = 're2_dist' + str(monomers) + case + '.dat'
fout = open(outfile, 'w')
for i in range(pdfRE2.size):
  minRE2 = i * dRE2
  maxRE2 = minRE2 + dRE2
  aveRE2 = minRE2 + dRE2 / 2.0
  outline = '%g %g %g %g\n' % (minRE2, maxRE2, aveRE2, pdfRE2[i])
  fout.write(outline)
fout.close()

sys.stdout.write('Output has been written to disk.\n')

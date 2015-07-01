"""This script loops over PNG files and adds a label to the
   bottom of each image using the convert utility of ImageMagick.
   Borders may be added as well.
"""

import sys
import os
import glob

# should a border to added ('yes' or 'no')
borders = 'yes'

# form a list of all PNG files
pngfiles = glob.glob('*.png')

# loop over files
for file in pngfiles:

  # determine values based on file name
  monomers = file[file.index('M') + 1:file.index('T')]
  timeval = '%.3f' % (int(file[file.index('T') + 1:file.index('I')]) / 1000000.0)
  RG = '%.1f' % (float(file[file.index('G') + 1:file.index('.p') - 1]))
  id = '%d' % (int(file[file.index('D') + 1:file.index('R')]))
  plane = file[file.index('.p') - 1:file.index('.p')]
  outfile = file

  # determine value of viewing plane
  if(plane == 'x'):
    plane = 'yz'
  elif(plane == 'y'):
    plane = 'zx'
  elif(plane == 'z'):
    plane = 'xy'
  else:
    sys.stdout.write('Problem with plane: ' + plane + '\n')

  # fixed labels
  fixed_label1 = '\' monomers per ring = ' + monomers + ', id = ' + id + ', viewing plane = ' + plane + \
                 '\n time = ' + timeval + ' million tau, RG = ' + RG + ' sigma\''

  # create and execute command
  cmd = 'convert ' + file + ' -background Khaki -font Courier-Bold -pointsize 14 ' + \
        'label:' + fixed_label1 + ' -gravity Center -append ' + file
  os.system(cmd)
  sys.stdout.write('Label added to ' + file + '\n')

  if(borders == 'yes'):
    # create and execute command
    cmd = 'convert -border 1x1 -bordercolor "#000000" ' + file + ' ' + outfile
    os.system(cmd)
    sys.stdout.write('Border added to ' + file + '\n')


"""
   This script loops over PNG files and adds a label to the
   bottom of each image using the convert utility of ImageMagick.
"""

import sys
import os
import glob

# form a list of all PNG files
pngfiles = glob.glob('*.png')

# loop over files
for file in pngfiles:

  # determine values based on file name
  monomers = file[file.index('M') + 1:file.index('T')]
  timeval = '%6.4f' % (int(file[file.index('T') + 1:file.index('R')]) / 1000000.0)
  outfile = file[0:file.index('R')] + 'label.png'

  # fixed labels
  fixed_label1 = '\' number of rings = 200, monomers per ring = ' + monomers
  fixed_label2 = '       time = '

  # create and execute command
  cmd = 'convert ' + file + ' -background Khaki -font Courier-Bold -pointsize 14 ' \
        'label:' + fixed_label1 + fixed_label2 + timeval + ' million tau\' ' \
        '-gravity Center -append ' + outfile
  os.system(cmd)
  sys.stdout.write('Label added to ' + file + '\n')

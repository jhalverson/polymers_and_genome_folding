"""
This script grabs all the PNG files in two directories and
make a single image out out of the corresponding pairs
(side-by-side).
"""

import sys
import os
import glob

pngfiles_left = glob.glob('../with_label/*.png')
pngfiles_rght = glob.glob('../with_label_rotated/*.png')

pngfiles_left.sort(lambda u, v: cmp(int(u[u.index('T') + 1:u.index('.png') - 5]), \
                                    int(v[v.index('T') + 1:v.index('.png') - 5])))
pngfiles_rght.sort(lambda u, v: cmp(int(u[u.index('T') + 1:u.index('.png') - 5]), \
                                    int(v[v.index('T') + 1:v.index('.png') - 5])))

for l,r in zip(pngfiles_left, pngfiles_rght):
  outfile = 'M1600' + l[l.index('T') + 1:l.index('.png') - 5] + 'side.png'
  cmd = 'convert ' + l + ' ' + r + ' +append ' + outfile
  os.system(cmd)

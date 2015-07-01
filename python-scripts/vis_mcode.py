"""This Python script sorts a set of still images and creates an animation
   of them using MEncoder."""

import os
import glob

blockfiles = glob.glob('*.png')
blockfiles = filter(lambda u: u.find('T') != -1, blockfiles)
blockfiles.sort(lambda u, v: cmp(int(u[u.index('T') + 1:u.index('.')]), \
                                 int(v[v.index('T') + 1:v.index('.')])))

s = ''
for file in blockfiles[:len(blockfiles) - 1]:
  s = s + file + ','
s = s + blockfiles[-1]

# high quality compression
cmd = 'mencoder "mf://' + s + '" -mf type=png:fps=24 -ovc lavc -lavcopts \
       vcodec=msmpeg4v2:vbitrate=8000:mbd=2:trell:vb_strategy=1:precmp=2:\
       cmp=2:subcmp=2 -o out.avi'

# same as Image Sequence in QuickTime Pro
cmd = 'mencoder "mf://' + s + '" -mf type=png:fps=24 -ovc copy -o out.avi'

failure = os.system(cmd)

#!/usr/bin/env python

"""This Python script gets files on a remote machine. It works by using
   ssh to create a file listing on the remote machine, downloading the
   listing and then getting only the files that it needs. These include
   new files and updated files. The script may be used to maintain single
   or multiple simulations. The core of the script may be replaced by
   rsync. Relative paths are used so the script should be run from the
   directory where it exists.

   Since cron provides a minimal environment it is necessary to set
   the PYTHONPATH manually.

   It may be advantageous to mirror the local and remote directory. That
   is, give them the same names and structure: ~/LOCAL/ptmp/jhalvers/foo
   versus the remote path of /ptmp/jhalvers/foo."""

import sys
import os
import glob

# format is [N, M_rings, M_linear, shear rates ... ]
N200M200M26 = ['200', '200', '26', '-6', '-5.5', '-5', '-4', '-3', '-2']
N400M200M26 = ['400', '200', '26', '-7', '-6.5', '-6', '-5.5', '-5', '-4', '-3', '-2']
cases = [N200M200M26, N400M200M26]
cases = [N400M200M26]

htmldir = '/people/thnfs/homes/halvers/public_html/rings/shear_contaminants'
keys = ['stress.dat.', 'velocity_profile.dat.', 'log.', 'restart.', 'in.', 'power6.sh', 'bluegene.sh']
scripts = ['viscosity', 'avestep', 'pij', 'vgrad', 'step']
imgs = ['viscosity', 'avestep', 'pij', 'vgrad', 'step', 'loop']

origdir = os.getcwd()
for c in cases:
  dr = 'N' + c[0] + 'M' + c[1] + 'M' + c[2]
  for r in range(len(c) - 3):
    drG = dr + 'G10' + c[r + 3]

    if(not os.path.isdir(drG)):
      sys.stdout.write('Making ' + drG + ' ...\n')
      os.mkdir(drG)
    os.chdir(drG)

    # create list of local files
    local = []
    for key in keys:
      local.extend(glob.glob(key + '*'))

    # create list of files on remote machine and copy it to local machine
    path = '/ptmp/jhalvers/sllod_' + dr + '/' + drG
    criteria = ' '.join([u + '*' for u in keys])
    os.system('ssh jhalvers@vip \'cd ' + path + '; rm -f listing.txt; ls -l ' + criteria + ' > listing.txt; exit\'')
    os.system('scp jhalvers@vip:' + path + '/listing.txt .')

    # download files if they are larger in size or not found on local machine
    f = open('listing.txt')
    data = f.readlines()
    f.close()
    os.system('rm -f listing.txt')
    new_files = False
    for d in data:
      if(any(key in d for key in keys)):
	fsz = d.split()[4]
        if(fsz.isdigit()):
          remote_file_size = int(fsz)
        else:
          remote_file_size = 10**15
          sys.stdout.write('\nWARNING: Did not find file size (' + fsz + ').\n')
	file = d.split()[-1]
        if(file in local):
          if(os.path.getsize(file) < remote_file_size):
            sys.stdout.write('Copying remote file (update): ' + file + '\n')
	    os.system('scp -p jhalvers@vip:' + path + '/' + file + ' .\n')
            new_files = True
        else:
          sys.stdout.write('Copying remote file (new): ' + file + '\n')
          os.system('scp -p jhalvers@vip:' + path + '/' + file + ' .\n')
          new_files = True

    # run analysis scripts and copy files to html directory
    if(new_files):
      for a in scripts:
        os.system('python ../figure_' + a + '.py')
      os.system('/bin/cp *.png ' + htmldir)

    os.chdir(origdir)
    sys.stdout.write('\n')

# create combined figure
os.system('python figure_combined.py')

# make HTML
s = """
<html>

<head>
<title>Shear viscosity of rings with linear</title>
</head>

<body>
11.1.2011
<p>
<h1>Shear viscosity of rings with linear</h1>

<table width="600"><tr><td>
Here we compute the viscosity as a function of shear rate for ring/linear blends.
</td></tr></table><p>
<img src="eta200.png" border="1"><p>
<img src="eta400.png" border="1"><p>
"""

for c in cases:
  dr = 'N' + c[0] + 'M' + c[1] + 'M' + c[2]
  s = s + '<h2>N = ' + c[0] + ', M_rings = ' + c[1] + ', M_linear = ' + c[2] + '</h2>\n'
  s = s + '<table><tr>\n'

  for img in imgs:
    s = s + '<tr>\n'
    for r in range(len(c) - 3):
      drG = dr + 'G10' + c[r + 3]
      s = s + '<td><img src="sllod_' + drG + '_' + img + '.png" border="1"></td>\n'
    s = s + '</tr>\n'

  s = s + '</table>\n'
s = s + '</body></html>'

f = open(htmldir + '/index.html', 'w')
f.write(s)
f.close()

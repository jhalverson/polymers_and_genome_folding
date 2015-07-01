#!/usr/bin/env python

"""Single simulation version.

   This Python script gets files on a remote machine. It works by using
   ssh to create a file listing on the remote machine, downloading the
   listing and then getting only the files that it needs. These include
   new files and updated files. The core of the script may be replaced by
   rsync. Relative paths are used so the script should be run from its
   own directory. The user will receive an email if a problem occurs in
   identifying file sizes or if a new directory for the HTML files is made.

   Since cron provides a minimal environment it is necessary to set
   the PYTHONPATH manually.

   Set htmldir to '' to turn off all html related features.

   It may be advantageous to mirror the local and remote directory. That
   is, give them the same names and structure: ~/LOCAL/ptmp/jhalvers/foo
   versus the remote path of /ptmp/jhalvers/foo. And the configuration
   files too."""

import sys
import os
import glob
import time

execfile('figure.common') # get job_name

htmldir = '/people/thnfs/homes/halvers/public_html/clc/' + job_name
remote_path = '/ptmp/jhalvers/explosion_ep_0.8T0.5G10-6'
keys = ['linear_', 'stress', 'velocity', 'log.', 'restart.', 'in.', 'power6.sh', 'bluegene.sh']
analysis_scripts = ['recon'] # assumed in ../scripts/
utility_scripts = ['step_loop', 'rg2re2 -p recon_linear_sllod_M200N200T0.5EP0.8G0.000001 -m 25 -x --Mrings 0 --Mlinear 200']
figure_scripts = ['rg2', 're2rg2', 'viscosity'] # assumed in ../scripts/
imgs = ['rg2', 're2rg2', 'viscosity', 'step', 'loop']

os.chdir('../configs')

####################################################################
##### IT SHOULD BE UNNECESSARY TO MAKE CHANGES BELOW THIS LINE #####
####################################################################

# create list of local files
local = []
for key in keys:
  local.extend(glob.glob(key + '*'))

# create list of files on remote machine and copy it to local machine
criteria = ' '.join([u + '*' for u in keys])
os.system('ssh jhalvers@vip \'cd ' + remote_path + '; rm -f listing.txt; ls -l ' + criteria + ' > listing.txt; exit\'')
os.system('scp jhalvers@vip:' + remote_path + '/listing.txt .')

# download files if they are larger in size or not found on local machine
f = open('listing.txt')
data = f.readlines()
f.close()
os.system('rm -f listing.txt')
trouble = ''
new_files = False
for d in data:
  if(any(key in d for key in keys)):
    file = d.split()[-1]
    fsz = d.split()[4]
    if(fsz.isdigit()):
      remote_file_size = int(fsz)
    else:
      remote_file_size = 10**15
      sys.stdout.write('\nWARNING: Did not find file size (' + fsz + ').\n')
      trouble += 'file size warning:' + file + ' ' + fsz + '\n'
    if(file in local):
      if(os.path.getsize(file) < remote_file_size):
	sys.stdout.write('Copying remote file (update): ' + file + '\n')
	os.system('scp -p jhalvers@vip:' + remote_path + '/' + file + ' .\n')
	new_files = True
    else:
      sys.stdout.write('Copying remote file (new): ' + file + '\n')
      os.system('scp -p jhalvers@vip:' + remote_path + '/' + file + ' .\n')
      new_files = True

if(trouble):
  cmd = 'echo "' + trouble + '" | mail -s "problem with ' + job_name + '" halverson@mpip-mainz.mpg.de'
  os.system(cmd)

if(new_files):
  # run analysis scripts, utilities and figure scripts
  for script in analysis_scripts:
    os.system('python ../scripts/' + script + '.py')
  for script in utility_scripts:
    os.system(script)
  for script in figure_scripts:
    os.system('python ../scripts/figure_' + script + '.py')

  if(htmldir != ''):
    # make directory
    if(not os.path.isdir(htmldir)):
      sys.stdout.write('Making ' + htmldir + ' ...\n')
      os.mkdir(htmldir)
      cmd = 'echo "' + htmldir + '" | mail -s "new html directory" halverson@mpip-mainz.mpg.de'
      os.system(cmd)

    # copy PNG files
    if(glob.glob('*.png') != []): os.system('/bin/cp *.png ' + htmldir)

    # make HTML
    s = '<html><head><title>' + job_name + '</title></head><body>\n'
    s += time.asctime() + '<p>\n'
    s += '<h1>' + job_name + '</h1>\n'
    s += '<table>\n'
    for img in imgs:
      s += '<tr><td><img src="' + job_name + '_' + img + '.png" border="1"></td></tr>\n'
    s += '</table></body></html>'
    f = open(htmldir + '/index.html', 'w')
    f.write(s)
    f.close()

    # update live jobs webpage
    f = open('/people/thnfs/homes/halvers/public_html/jobs.html')
    data = f.readlines()
    f.close()
    if(all([job_name not in d for d in data])):
      data.insert(2, '<a href="' + htmldir[htmldir.index('public_html/') + 12:] + '/index.html">' + job_name + '</a><br>\n')
      f = open('/people/thnfs/homes/halvers/public_html/jobs.html', 'w')
      f.writelines(data)
      f.close()

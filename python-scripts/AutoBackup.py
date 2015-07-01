#!/usr/bin/env python

"""This Python script backs-up files from simulation jobs. It works by
   going to the job directory, creating an inventory of all the files,
   compares this to the previous inventory and adds the new files that
   are more than 24 hrs old to a tar file. The step bounds and job step
   are saved in the filename along with the date.

   Possible improvement would be to walk the directory tree and determine
   automatically if the directory should be backed-up.

   A flag could be added which would allow one to store files that have
   been modified since archived. Currently only the original is stored.
"""

import sys
import os
import time
import shutil
import tarfile

p = '/ptmp/jhalvers/'
r = '/r/j/jhalvers/'
email = 'jhalvers@rzg.mpg.de'
step_determination_ignore = ['log.', 'stress.', 'velocity_profile.']
archive_ignore = ['listing.txt', 'archived.files']

# format of dr is [job directory on /ptmp, tar file directory on /r]
dr = []
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-7',   r + 'new_shear'])
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-6.5', r + 'new_shear'])
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-6',   r + 'new_shear'])
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-5.5', r + 'new_shear'])
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-5',   r + 'new_shear'])
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-4',   r + 'new_shear'])
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-3',   r + 'new_shear'])
dr.append([p + 'sllod_N400M200M26/N400M200M26G10-2',   r + 'new_shear'])

#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-7',   r + 'new_shear'])
#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-6.5', r + 'new_shear'])
#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-6',   r + 'new_shear'])
#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-5.5', r + 'new_shear'])
#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-5',   r + 'new_shear'])
#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-4',   r + 'new_shear'])
#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-3',   r + 'new_shear'])
#dr.append([p + 'sllod_N400M113M113/N400M113M113G10-2',   r + 'new_shear'])

dr.append([p + 'N200_M169_M57_INT32_3',   r + 'new_contam'])
dr.append(['/u/jhalvers/BlueGene/N400M200M3', r + 'new_contam'])

dr.append([p + 'sllod_N200M200M6/N200M200M6G10-6',   r + 'new_shear'])
dr.append([p + 'sllod_N200M200M6/N200M200M6G10-5.5', r + 'new_shear'])
dr.append([p + 'sllod_N200M200M6/N200M200M6G10-5',   r + 'new_shear'])
dr.append([p + 'sllod_N200M200M6/N200M200M6G10-4.5', r + 'new_shear'])
dr.append([p + 'sllod_N200M200M6/N200M200M6G10-4',   r + 'new_shear'])
dr.append([p + 'sllod_N200M200M6/N200M200M6G10-3',   r + 'new_shear'])
dr.append([p + 'sllod_N200M200M6/N200M200M6G10-2',   r + 'new_shear'])
dr.append([p + 'sllod_N200M200M6/N200M200M6_init',   r + 'new_shear'])

dr.append([p + 'sllod_N200M200M26/N200M200M26G10-5.75', r + 'new_shear'])
dr.append([p + 'sllod_M200N800G10-6.5', r + 'new_shear'])

dr.append([p + 'sllod_N200M100M100/N200M100M100G10-6',   r + 'new_shear'])
dr.append([p + 'sllod_N200M100M100/N200M100M100G10-5.5', r + 'new_shear'])
dr.append([p + 'sllod_N200M100M100/N200M100M100G10-5',   r + 'new_shear'])
dr.append([p + 'sllod_N200M100M100/N200M100M100G10-4',   r + 'new_shear'])
dr.append([p + 'sllod_N200M100M100/N200M100M100G10-3',   r + 'new_shear'])
dr.append([p + 'sllod_N200M100M100/N200M100M100G10-2',   r + 'new_shear'])
dr.append([p + 'sllod_N200M100M100/N200M100M100_init',   r + 'new_shear'])

dr.append([p + 'explosion_ep_0.8T0.6G10-6', r + 'clc'])
dr.append([p + 'explosion_ep_0.8T0.5G10-6', r + 'clc'])
dr.append([p + 'explosion_ep_0.8T0.6G0', r + 'clc'])
dr.append([p + 'explosion_ep_0.8T0.5G0', r + 'clc'])

####################################################################
##### IT SHOULD BE UNNECESSARY TO MAKE CHANGES BELOW THIS LINE #####
####################################################################

char72 = '123456789012345678901234567890123456789012345678901234567890123456789012\n\n'
s_per_day = 24 * 60 * 60

msg = {}
msg[0] = 'Directory on /ptmp does not exist: '
msg[1] = 'Directory on /r does not exist: '
msg[2] = 'Starting backup of '
msg[3] = 'AutoBackup: File already exists'
msg[4] = 'AutoBackup: Directory does not exist'
msg[5] = 'AutoBackup: New archive'
msg[6] = 'AutoBackup: Summary'

def min_max_step(files):
  """Find smallest and largest step."""
  step = []
  for f in files:
    if(all([u not in f for u in step_determination_ignore]) and '.' in f):
      suffix = f[f.rindex('.') + 1:]
      if(suffix.isdigit()):
        step.append(int(suffix))
  if(len(step) > 0):
    mn = min(step)
    mn = str(mn / 10**5) + 'E5' if mn >= 10**5 else str(mn)
    mx = max(step)
    mx = str(mx / 10**5) + 'E5' if mx >= 10**5 else str(mx)
    return (mn, mx)
  else:
    return ('_', '_')

def send_email(body, subject, address):
   """Send email message."""
   os.system('echo "' + char72 + body + '" | mail -s "' + subject + '" ' + address)

def create_archive(files, m, pdir, rdir):
  """"Create archive and move to backup file system."""
  now = time.strftime('%d%b%Y', time.localtime())
  jobstep = '_'
  if(os.path.isfile('jobstep.dat')):
    f = open('jobstep.dat')
    jobstep = '_J' + f.readline().strip() + '_'
    f.close()
  flnm = os.path.split(pdir)[1] + '_' + m[0]  + '_' + m[1] + jobstep + now + '.tar'
  tar = tarfile.open(flnm, 'w:')
  for f in files:
    tar.add(f)
  tar.close()
  if(os.path.isfile(os.path.join(rdir, flnm))):
    send_email(f + ': File not backed-up.', msg[3], email)
    return flnm + ': --'
  else:
    file_size = str(os.path.getsize(flnm))
    shutil.copy2(flnm, os.path.join(rdir, flnm))
    os.remove(flnm)
    return flnm + ': Number of files = ' + str(len(files)) + ', TAR size = ' + file_size + ' bytes\n'

summary = ''
for pdir, rdir in dr:
  sys.stdout.write('Working on ' + pdir + ' ...\n')
  if(not os.path.isdir(pdir)):
    send_email(msg[0] + pdir, msg[4], email)
  elif(not os.path.isdir(rdir)):
    send_email(msg[1] + rdir, msg[4], email)
  else:
    os.chdir(pdir)
    all_files = os.listdir(os.curdir)
    if(os.path.isfile('archived.files')):
      f = open('archived.files'); archived_files = eval(f.readline()); f.close()
      files = filter(lambda u: time.time() - os.path.getmtime(u) > s_per_day and \
                               os.path.isfile(u) and \
                               u not in archived_files and \
                               u not in archive_ignore, all_files)
      f = open('archived.files', 'w'); f.write(repr(archived_files + files)); f.close()
    else:
      files = filter(lambda u: time.time() - os.path.getmtime(u) > s_per_day and \
                               os.path.isfile(u) and \
                               u not in archive_ignore, all_files)
      if(len(files) > 0):
        f = open('archived.files', 'w'); f.write(repr(files)); f.close()
        send_email(msg[2] + pdir, msg[5], email)
    if(len(files) > 0): summary += create_archive(files, min_max_step(files), pdir, rdir)

# send summary email
send_email(summary, msg[6], email)

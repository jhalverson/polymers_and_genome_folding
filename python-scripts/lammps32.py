"""This Python script changes file names and time step values in
   the file so that files from the post-INT32 job are continuous
   with the pre-INT32 job. The user should manually check that
   the last file from pre-INT32 is identical to the 0 file of
   post-INT32. And then remove the 0 file.

   The idea is to run the script from the same directory where
   the files are that need their step values change. Then move
   the corrected files to the desired directory.

   Note that the file list is reverse sorted. This is to ensure
   that when the file is corrected it does not overwrite a file
   that has yet to be corrected which can happen if last_timestep
   is smaller than the overall time of the file set and the list
   is not reverse sorted.
"""

import sys
import os
import glob

last_timestep = 2050000000 # last timestep of pre-INT32 files
postINT32_prefix = 'rings200_3' # prefix of files in current directory to be corrected
preINT32_prefix = '../configs/200_3' # path and prefix to files with correct step values

files = glob.glob(postINT32_prefix + '.*')
files = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), files)
files.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])), reverse=True)
if(files.count(postINT32_prefix + '.0')):
  print 'remove t=0 file from post-INT32 set and check that it is the same the first file in the pre set'
  sys.exit()

for file in files:
  f = open(file)
  data = f.readlines()
  f.close()

  # correct the time step
  corrected_timestep = int(data[1]) + last_timestep
  data[1] = str(corrected_timestep) + '\n'

  # write corrected file
  outfile = postINT32_prefix + '.' + str(corrected_timestep)
  f = open(outfile, 'w')
  f.writelines(data)
  f.close()

  # remove uncorrected file
  ec = os.system('rm -f ' + file)

# consistency check by file name
files = glob.glob(preINT32_prefix + '.*') + glob.glob(postINT32_prefix + '.*')
files = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), files)
files.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))

exts = []
for i in range(len(files)):
  file_extension = int(files[i][files[i].rindex('.') + 1:])
  exts.append(file_extension)

for i in range(len(files)):
  if(exts.count(i * 100000)) == 0:
    print 'missing file with extension:', i * 100000

for i in range(len(files)):
  file_extension = int(files[i][files[i].rindex('.') + 1:])
  if(file_extension != i * 100000):
    print 'problem', file_extension
    sys.exit(1)

print 'Files found to be consistent'
print files[0]
print files[1]
print files[-1]

"""This script computes g3 for each chain as a function of time. It is
   not possible to average over time origins. In output file chain indices
   range from 1 to chains."""

# use sys for printing compatibility between Python 2 & 3
import sys
import os
import glob
import pypolymer

# M_rings + M_linear
chains = 206

# generation of file names: (1) glob, (2) loop or (3) manual
blockfiles_method = 1

if(blockfiles_method == 1):

  prefix = '../../rings400_6_slow_insert.'
  existingfiles = glob.glob(prefix + '*') 
  existingfiles = filter(lambda u: u[u.rindex('.') + 1:].isdigit(), existingfiles)
  starting_timestep = 100000000
  existingfiles = filter(lambda u: int(u[u.rindex('.') + 1:]) >= starting_timestep, existingfiles)
  existingfiles.sort(lambda u, v: cmp(int(u[u.rindex('.') + 1:]), int(v[v.rindex('.') + 1:])))
  blockfiles = existingfiles[::]                                                              
  if 1:                                                                                       
    maxfiles = 20                                                                            
    if(len(blockfiles) > maxfiles):                                                           
      blockfiles = blockfiles[::int(len(blockfiles) / float(maxfiles))]   

elif(blockfiles_method == 2):

  istr = 10000000
  iend = 22500000
  incr = 500000

  blockfiles = []
  for i in range(istr, iend + incr, incr):
    blockfiles.append('RG_N200M1600.0' + str(i))

elif(blockfiles_method == 3):

  blockfiles = ['RG_N200M1600.0361200000', 'RG_N200M1600.0361300000']
 
else:

  sys.stdout.write('ERROR: blockfiles_method\n')
  sys.exit(1)

if(len(blockfiles) == 0):
  sys.stderr.write('ERROR: blockfiles is empty\n')
  sys.exit()                                      

sys.stdout.write('\n\n')
sys.stdout.write('=================================================\n')
sys.stdout.write('                     g3_history                  \n')
sys.stdout.write('=================================================\n')
sys.stdout.write('\n\nTotal number of files: ' + str(len(blockfiles)) + '\n')
if(len(blockfiles) > 2):                                                                  
  sys.stdout.write('  blockfiles[0] = ' + blockfiles[0] + '\n')                           
  sys.stdout.write('  blockfiles[1] = ' + blockfiles[1] + '\n')                           
  sys.stdout.write('  blockfiles[2] = ' + blockfiles[2] + '\n')                           
  sys.stdout.write('  ...' + '\n')                                                        
  sys.stdout.write('  blockfiles[-3] = ' + blockfiles[-3] + '\n')                         
  sys.stdout.write('  blockfiles[-2] = ' + blockfiles[-2] + '\n')                         
  sys.stdout.write('  blockfiles[-1] = ' + blockfiles[-1] + '\n\n')  

# which chains should be written to file: (1) all, (2) every nth or (3) only a few
chains2file = 1

if(chains2file == 1):
  # all chain data will be written
  indices = range(chains)
elif(chains2file == 2):
  # every nth chain
  nth = 20
  indices = range(0, chains, nth)
elif(chains2file == 3):
  # only these chains will have their data written
  indices = range(200, chains)
else:
  sys.stdout.write('Value of chains2file is not valid.\n\n')
  sys.exit(1)

# initialize the list
time_g3 = []

# initialize COM list at t = 0
com_zero = []

# loop over files
for file in blockfiles:

  # initialize the sublist
  time_g3_sub = []

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
    time_zero = timestep
  if(total_particles1st != total_particles):                  
    sys.stderr.write('Change in total_particles (' + str(total_particles) + ').\n')
    sys.exit()                                                                     

  sys.stdout.write('Working on '+ file + '\n')

  # store the time value (time steps)
  time_g3_sub.append(timestep)

  # use first file to get initial positions
  if(file == blockfiles[0]):
    xcm_sys_zero = 0.0; ycm_sys_zero = 0.0; zcm_sys_zero = 0.0
    for i in range(chains):
      xcm = 0.0; ycm = 0.0; zcm = 0.0
      for j in range(monomers):
        p = i * monomers + j
        xcm = xcm + x[p]
        ycm = ycm + y[p]
        zcm = zcm + z[p]
        xcm_sys_zero = xcm_sys_zero + x[p]
	ycm_sys_zero = ycm_sys_zero + y[p]
	zcm_sys_zero = zcm_sys_zero + z[p]
      xcm = xcm / monomers
      ycm = ycm / monomers
      zcm = zcm / monomers

      # store initial center-of-mass position
      com_zero.append((xcm, ycm, zcm))

    xcm_sys_zero = xcm_sys_zero / total_particles
    ycm_sys_zero = ycm_sys_zero / total_particles
    zcm_sys_zero = zcm_sys_zero / total_particles
    sys.stdout.write('r_COM(t=' + str(timestep) + '): %8.2f %8.2f %8.2f\n' % \
                    (xcm_sys_zero, ycm_sys_zero, zcm_sys_zero))

  xcm_sys = 0.0; ycm_sys = 0.0; zcm_sys = 0.0
  for i in range(chains):
    xcm = 0.0; ycm = 0.0; zcm = 0.0
    for j in range(monomers):
      p = i * monomers + j
      xcm = xcm + x[p]
      ycm = ycm + y[p]
      zcm = zcm + z[p]
      xcm_sys = xcm_sys + x[p]
      ycm_sys = ycm_sys + y[p]
      zcm_sys = zcm_sys + z[p]
    xcm = xcm / monomers
    ycm = ycm / monomers
    zcm = zcm / monomers

    # compute and store mean squared displacement
    msd = (xcm - com_zero[i][0])**2 + (ycm - com_zero[i][1])**2 + (zcm - com_zero[i][2])**2
    time_g3_sub.append(msd)

  xcm_sys = xcm_sys / total_particles
  ycm_sys = ycm_sys / total_particles
  zcm_sys = zcm_sys / total_particles
  if(abs(xcm_sys-xcm_sys_zero) > 0.01 or \
     abs(ycm_sys-ycm_sys_zero) > 0.01 or \
     abs(zcm_sys-zcm_sys_zero) > 0.01):
    sys.stdout.write('WARNING: r_COM(t=%d steps): %8.2f %8.2f %8.2f\n' % \
    (timestep,xcm_sys-xcm_sys_zero,ycm_sys-ycm_sys_zero,zcm_sys-zcm_sys_zero))

  # append sublist to list
  time_g3.append(time_g3_sub)

# write out the data
outfile = 'g3_history_' + str(M_rings) + '_' + str(M_linear) + '_N' + str(monomers) + '.dat'
fout = open(outfile, 'w')

# create comment line
s = ''
for index in indices:
  s = s + 'i=%-3d     ' % (index + 1)

# write comment line
fout.write('  # timestep   '  + s + '\n')

# write out data
for i in range(len(blockfiles)):
  s = ''
  for index in indices:
    s = s + '%9.2f ' % (time_g3[i][index + 1])
  time_fmt = '%10d ' % (time_g3[i][0] - time_zero)
  fout.write(time_fmt + s + '\n')
fout.close()
sys.stdout.write(outfile + ' has been written to disk (g3 vs. time data file).\n')

"""This script determines which configuration files should
   be chosen to smoothly show the ppa analysis which has
   large changes at first and then slow changes. The idea
   is that the difference in sucessive configuration files
   should increase with the frame number."""

# frames per second
frame_rate = 24

# total length of animation in seconds
animation_time = 24

# number of frames to be generated
total_frames = frame_rate * animation_time

# time step at first configuration
cfg_start = 0

# time step at last configuration
cfg_last = 680000

# a function that determines the difference between successive configs at each frame
def f(frame_index):
  return frame_index**1.2395

print 'frame  time_step'
print 1, cfg_start

# init
cmd = []

prev = cfg_start
for i in range(2, total_frames + 1):
  nxt = f(i)
  fi = prev + nxt
  prev = fi

  # throw away precision if necessary
  if(10000 > int(fi)):
    a = int(fi)
    s = 'scp jhalvers@vip:/ptmp/jhalvers/ppa_rings/lammps_ppa800_bg_early_256/ppa800a.' + str(a) + ' .\n'
    cmd.append(s)
  elif(10000 < int(fi) < 385800):
    a = int(fi) / 100
    a = a * 100
    s = 'scp jhalvers@vip:/ptmp/jhalvers/ppa_rings/lammps_ppa800_bg_long_256/ppa800a.' + str(a) + ' .\n'
    cmd.append(s)
  else:
    a = int(fi) / 1000
    a = a * 1000
    s = 'scp jhalvers@vip:/ptmp/jhalvers/ppa_rings/lammps_ppa800_bg_long_256_contd/ppa800a.' + str(a) + ' .\n'
    cmd.append(s)

  print i, a, int(nxt)

outfile = open('scp_commands.sh', 'w')
outfile.write('scp jhalvers@vip:/ptmp/jhalvers/ppa_rings/lammps_ppa800_bg_early_256/ppa800a.0 .\n')
for j in range(len(cmd)):
  outfile.write(cmd[j])
outfile.close()

"""This Python script is useful for generating groups
   in LAMMPS to do a PPA."""

import sys
import random

M_rings = 200
M_linear = 0
N = 1600

# create list of fixed particles
N2 = N / 2

def get_rnd():
  num = int(N2 * random.random()) + 1
  if(num > N2 or num < 1): sys.exit()
  return num

# rings
ids = []
for i in range(M_rings):
  rnum = get_rnd()
  tmp = i * N + rnum
  ids.append(tmp)
  ids.append(tmp + N2)
  chain_id1 = (tmp - 1) / N
  chain_id2 = (tmp + N2 - 1) / N
  if(chain_id1 != chain_id2): print i*N, i*N+N, rnum, tmp, tmp+N2, chain_id1, chain_id2

# linear
for i in range(0, M_linear):
  ids.append(M_rings * N + 1 + i * N)
  ids.append(M_rings * N + 1 + i * N + N - 1)

# reverse the list
ids = ids[::-1]

# shoud be 400 ids in list
for i in range(25):
  s = 'group           cntMove' + str(i+1) + ' id '
  for j in range(16):
    s = s + str(ids.pop()) + ' '
  print s

# linear
if(M_linear != 0):
  i = 25
  s = 'group           cntMove' + str(i+1) + ' id '
  for j in range(12):
    s = s + str(ids.pop()) + ' '
  print s

import numpy as np

# store list
lst = dir(np)

# how many columns
c = 2

# use natural lengths (-1) or impose a maximum value (#)
length = -1

if(length == -1):
  max = 0
  for i in range(0, len(lst)):
    if(len(lst[i]) > max): max = len(lst[i])
else:
  max = length

# build string of spaces
emp = ''
for i in range(max):
 emp = emp + ' '

# number of rows
rows = len(lst) / c
extra = len(lst) % c

def add_spaces(x):
  if(length == -1):
    return x + emp[:max - len(x) + 1]
  else:
    if(len(x) < max):
      x = x + emp[:max - len(x) + 1]
      return x[:max - 2] + '|'
    else:
      return x[:max - 3] + '.|'

for i in range(rows):
  s = ''
  for j in range(c):
    s = s  + add_spaces(lst[(i * c + j)])
  print s

if(extra != 0):
  s = ''
  for j in range(extra):
    s = s  + add_spaces(lst[(i * c + j)])
  print s

"""This Python script launches several jobs with the same          
   starting configuration but different shear rates."""            

import os

N = '200'
M_rings = '200'
M_linear = '26'
shear_rates = ['-5.5', '-5', '-4', '-3', '-2']

srDict = {}
srDict['-7'] =   '0.0000001'
srDict['-6.5'] = '0.000000316227766'
srDict['-6'] =   '0.000001'         
srDict['-5.5'] = '0.00000316227766'
srDict['-5'] =   '0.00001'
srDict['-4.5'] = '0.0000316227766'
srDict['-4'] =   '0.0001'
srDict['-3.5'] = '0.000316227766'
srDict['-3'] =   '0.001'
srDict['-2'] =   '0.01'

origdir = os.getcwd()
for r in shear_rates:
  dr = 'N' + N + 'M' + M_rings + 'M' + M_linear
  drG = dr + 'G10' + r
  if(not os.path.exists(drG)):
    print 'Making ' + drG + ' ...'
    os.mkdir(drG)
    os.system('cp ' + dr + '_init/in.sllod ' + drG)
    os.system('cp ' + dr + '_init/jobstep.dat ' + drG)
    os.system('cp ' + dr + '_init/bluegene.sh ' + drG)

    # modify in.sllod
    f = open(drG + '/in.sllod')
    data = f.readlines()
    f.close()
    data[2] = 'variable        shear_rate string ' + srDict[r] + '\n'
    f = open(drG + '/in.sllod', 'w')
    f.writelines(data)
    f.close()

    # modify bluegene.sh
    f = open(drG + '/bluegene.sh')
    data = f.readlines()
    f.close()
    data[1] = '# @ job_name = ' + drG  + '\n'
    f = open(drG + '/bluegene.sh', 'w')
    f.writelines(data)
    f.close()

    # submit job
    os.chdir(drG)
    os.system('llsubmit bluegene.sh')
    os.chdir(origdir)

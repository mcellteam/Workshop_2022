import os,sys
import subprocess as sp
import bionetgen
import numpy as np
import matplotlib.pyplot as plt

'''This script runs mcell & bionetgen from the same bngl model & plot the results together
GCG
02.10.22
'''
#
python_path ='/Applications/Blender-2.93-CellBlender/blender.app/Contents/Resources/2.93/python/bin/python3.9'
run_file = 'model.py'
tdir = './'
# # # # #run mcell
Proc = sp.call([python_path, run_file],cwd = tdir)
if Proc != 0:
    print('MCell did not run')
else:
    print('MCell sim done')

# for Windows:
# replace C:\\cmw2021\\ with the actual path where you unpacked CellBlender
# python_path = 'C:\\cmw2021\\Blender-2.93-CellBlender\\2.93\\python\\bin\\python3.9.exe'


#Load mcell output
mc_snare_syn = np.genfromtxt('./react_data/seed_00001/SNARE_sync.dat',
                      dtype=float,
                      delimiter=' ')
#
mc_snare_asyn = np.genfromtxt('./react_data/seed_00001/SNARE_async.dat',
                      dtype=float,#
                      delimiter=' ')
mc_vrel = np.genfromtxt('./react_data/seed_00001/V_release.dat',
                      dtype=float,#
                      delimiter=' ')
rel_glu = np.genfromtxt('./glu_rel.dat',
                      dtype=float,#
                      delimiter=' ')
print(len(rel_glu), np.shape(rel_glu))
fig, ax = plt.subplots()
fig.subplots_adjust(right=0.9, left = 0.15, bottom =0.15, top = 0.95)
#
#ax.plot(bng[:,0],bng[:,2],'r',label = 'BNGL ODE SNARE_sync')
#ax.plot(bng[:,0],bng[:,4],'g',label = 'BNGL SNARE_async')
#ax.plot(bng[:,0],bng[:,6],'b',label = 'BNGL ODE V release')
ax.plot(mc_snare_syn[:,0],mc_snare_syn[:,1],'r',linestyle = '--',label = 'MCell SNARE_sync')
ax.plot(mc_snare_asyn[:,0],mc_snare_asyn[:,1],'g',linestyle = '--',label = 'MCell SNARE_async')
ax.plot(mc_vrel[1:,0],np.diff(mc_vrel[:,1]),'b',linestyle = '--',label = 'MCell V release')
#ax.plot(mc_vrel[:,0],mc_vrel[:,1],'b',linestyle = '--',label = 'MCell V release')
ax.plot(rel_glu[:,0],np.ones(len(rel_glu[:,0])),'. m', ms = '10')

plt.legend(loc='upper left')
#plt.ylim(0,8000)
#plt.xlim(0,3.5e-4)
#
plt.xlabel('Time (sec)')
plt.ylabel('# Molecules')
#plt.savefig('snare.png',dpi=600)
plt.show()

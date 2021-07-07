#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loadeddata#
parameter=np.loadtxt('../../Outfile/Parameter.dat')

#constants#
enetrans = 5.59429*10**(-55)
timetrans = 2.03017*10**(5)
lengthtrans = 1.4766839
rhotrans = 1.6190*10**(-18)

#load#
dx2 = parameter[1, 1]
length2 = parameter[1, 0]

#grid#
grid2 = int(length2/dx2 + 10)
number2 = grid2 + 2

#loaddata#
scag=np.loadtxt('../../Outfile/Flame/Star_WENO_ScaG_1.dat',comments='"')
scag2=np.loadtxt('../../Outfile/Flame/Star_WENO_ScaG2_1.dat',comments='"')

# get grid #
nmcount = int(np.size(scag[:,0])/grid2)

#array#
nmtime = np.zeros(nmcount)

# load time array #
for i in range(0, nmcount):
    nmtime[i] = np.genfromtxt('../../Outfile/Hydro/Star_WENO_Density_NM_1.dat', skip_header = i * (grid2 + 2), max_rows = 1, usecols = (2))

#plot#
for i in range(0, nmcount):
   plt.plot(scag[i*grid2 + 5:(i+1)*grid2 - 5,0], scag[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Level-Sset (CODE)')
plt.title('Flame Profile')
plt.grid(True)
plt.savefig('flame.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(scag2[i * grid2 + 5:(i + 1) * grid2 - 5, 0], scag2[i * grid2 + 5:(i + 1) * grid2 - 5, 1], label='%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Level-Sset (CODE)')
plt.title('Deton Profile')
plt.grid(True)
plt.savefig('deton.png')
plt.clf()

#prevent run-time error#
plt.close('all')
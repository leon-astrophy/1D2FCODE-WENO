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
rho2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Density_NM_1.dat',comments='"')
eps2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Epsilon_NM_1.dat',comments='"')
pre2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Pressure_NM_1.dat',comments='"')
vel2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Velocity_NM_1.dat',comments='"')

# get grid #
nmcount = int(np.size(rho2[:,0])/grid2)

#array#
nmtime = np.zeros(nmcount)

# load time array #
for i in range(0, nmcount):
    nmtime[i] = np.genfromtxt('../../Outfile/Hydro/Star_WENO_Density_NM_1.dat', skip_header = i * (grid2 + 2), max_rows = 1, usecols = (2))

#plot#
for i in range(0, nmcount):
   plt.plot(rho2[i*grid2 + 5:(i+1)*grid2 - 5,0], rho2[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Density (CODE)')
plt.title('Density Profile')
plt.grid(True)
plt.savefig('rho2.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(eps2[i*grid2 + 5:(i+1)*grid2 - 5,0], eps2[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Epsilon (CODE)')
plt.title('Epsilon Profile')
plt.grid(True)
plt.savefig('eps2.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(pre2[i*grid2 + 5:(i+1)*grid2 - 5,0], pre2[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Pressure (CODE)')
plt.title('Pressure Profile')
plt.grid(True)
plt.savefig('pre2.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(vel2[i*grid2 + 5:(i+1)*grid2 - 5,0], vel2[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Velocity (CODE)')
plt.title('Velocity Profile')
plt.grid(True)
plt.savefig('vel2.png')
plt.clf()

#prevent run-time error#
plt.close('all')
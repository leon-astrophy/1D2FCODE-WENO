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
dx1 = parameter[0, 1]
length1 = parameter[0, 0]

#grid#
grid1 = int(length1/dx1 + 10)
number1 = grid1 + 2

#loaddata#
rho1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Density_DM_1.dat',comments='"')
eps1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Epsilon_DM_1.dat',comments='"')
pre1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Pressure_DM_1.dat',comments='"')
vel1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Velocity_DM_1.dat',comments='"')

# get grid #
dmcount = int(np.size(rho1[:,0])/grid1)

#array#
dmtime = np.zeros(dmcount)

# load time array #
for i in range(0, dmcount):
    dmtime[i] = np.genfromtxt('../../Outfile/Hydro/Star_WENO_Density_DM_1.dat', skip_header = i * (grid1 + 2), max_rows = 1, usecols = (2))

#plot#
for i in range(0, dmcount):
   plt.plot(rho1[i*grid1 + 5:(i+1)*grid1 - 5,0], rho1[i*grid1 + 5:(i+1)*grid1 - 5,1], label = '%.3f' % dmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Density (CODE)')
plt.title('Density Profile')
plt.grid(True)
plt.savefig('rho1.png')
plt.clf()

#plot#
for i in range(0, dmcount):
   plt.plot(eps1[i*grid1 + 5:(i+1)*grid1 - 5,0], eps1[i*grid1 + 5:(i+1)*grid1 - 5,1], label = '%.3f' % dmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Epsilon (CODE)')
plt.title('Epsilon Profile')
plt.grid(True)
plt.savefig('eps1.png')
plt.clf()

#plot#
for i in range(0, dmcount):
   plt.plot(pre1[i*grid1 + 5:(i+1)*grid1 - 5,0], pre1[i*grid1 + 5:(i+1)*grid1 - 5,1], label = '%.3f' % dmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Pressure (CODE)')
plt.title('Pressure Profile')
plt.grid(True)
plt.savefig('pre1.png')
plt.clf()

#plot#
for i in range(0, dmcount):
   plt.plot(vel1[i*grid1 + 5:(i+1)*grid1 - 5,0], vel1[i*grid1 + 5:(i+1)*grid1 - 5,1], label = '%.3f' % dmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Velocity (CODE)')
plt.title('Velocity Profile')
plt.grid(True)
plt.savefig('vel1.png')
plt.clf()


#prevent run-time error#
plt.close('all')
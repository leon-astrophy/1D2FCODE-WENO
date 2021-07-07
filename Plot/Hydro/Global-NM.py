#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
mass2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Mass_NM_1.dat',comments='"')
energy2=np.loadtxt('../../Outfile/Hydro/Star_WENO_Energy_NM_1.dat',comments='"')
intenergy2=np.loadtxt('../../Outfile/Hydro/Star_WENO_IntEnergy_NM_1.dat',comments='"')
mechenergy2=np.loadtxt('../../Outfile/Hydro/Star_WENO_MechEnergy_NM_1.dat',comments='"')
centralrho2=np.loadtxt('../../Outfile/Hydro/Star_WENO_CentralDensity_NM_1.dat',comments='"')

#plot#
plt.plot(centralrho2[:,0], centralrho2[:,1], label='Central Density')
plt.plot(centralrho2[:,0], centralrho2[:,2], label='Maximum Density')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Density (CODE)')
plt.title('Density Versus Time')
plt.grid(True)
plt.savefig('rhoc2.png')
plt.clf()

#plot#
plt.plot(mass2[:,0], mass2[:,1])
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Mass (CODE)')
plt.title('Mass Versus Time')
plt.grid(True)
plt.savefig('mass2.png')
plt.clf()

#plot#
plt.plot(energy2[:,0], energy2[:,1], label='Total Energy')
plt.plot(intenergy2[:,0], intenergy2[:,1], label='Internal Energy')
plt.plot(mechenergy2[:,0], mechenergy2[:,1], label='Kinetic Energy')
plt.plot(mechenergy2[:,0], mechenergy2[:,2], label='Gravitational Energy')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy (CODE)')
plt.title('Energy Versus Time')
plt.grid(True)
plt.savefig('energy2.png')
plt.clf()

#plot#
plt.plot(energy2[:,0], energy2[:,2])
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy (CODE)')
plt.title('Energy Input Versus Time')
plt.grid(True)
plt.savefig('energyinput.png')
plt.clf()

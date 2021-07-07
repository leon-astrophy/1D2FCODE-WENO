#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
mass1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Mass_DM_1.dat',comments='"')
energy1=np.loadtxt('../../Outfile/Hydro/Star_WENO_Energy_DM_1.dat',comments='"')
intenergy1=np.loadtxt('../../Outfile/Hydro/Star_WENO_IntEnergy_DM_1.dat',comments='"')
mechenergy1=np.loadtxt('../../Outfile/Hydro/Star_WENO_MechEnergy_DM_1.dat',comments='"')
centralrho1=np.loadtxt('../../Outfile/Hydro/Star_WENO_CentralDensity_DM_1.dat',comments='"')

#plot#
plt.plot(centralrho1[:,0], centralrho1[:,1], label='Central Density')
plt.plot(centralrho1[:,0], centralrho1[:,2], label='Maximum Density')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Density (CODE)')
plt.title('Density Versus Time')
plt.grid(True)
plt.savefig('rhoc1.png')
plt.clf()

#plot#
plt.plot(mass1[:,0], mass1[:,1])
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Mass (CODE)')
plt.title('Mass Versus Time')
plt.grid(True)
plt.savefig('mass1.png')
plt.clf()

#plot#
plt.plot(energy1[:,0], energy1[:,1], label='Total Energy')
plt.plot(intenergy1[:,0], intenergy1[:,1], label='Internal Energy')
plt.plot(mechenergy1[:,0], mechenergy1[:,1], label='Kinetic Energy')
plt.plot(mechenergy1[:,0], mechenergy1[:,2], label='Gravitational Energy')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy (CODE)')
plt.title('Energy Versus Time')
plt.grid(True)
plt.savefig('energy1.png')
plt.clf()

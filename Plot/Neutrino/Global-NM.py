#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
pair=np.loadtxt('../../Outfile/Neutrino/Star_WENO_PairNeutrinoSpectra_1.dat',comments='"')
plas=np.loadtxt('../../Outfile/Neutrino/Star_WENO_PlasNeutrinoSpectra_1.dat',comments='"')
thermo=np.loadtxt('../../Outfile/Neutrino/Star_WENO_ThermoNeutrinoEnergyLoss_1.dat',comments='"')

#plot#
plt.plot(pair[:,0], pair[:,1], label='1 MeV')
plt.plot(pair[:,0], pair[:,2], label='2 MeV')
plt.plot(pair[:,0], pair[:,3], label='3 MeV')
plt.plot(pair[:,0], pair[:,4], label='4 MeV')
plt.plot(pair[:,0], pair[:,5], label='5 MeV')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Production Rate (CODE)')
plt.title('Pair Neutrino Spectra')
plt.grid(True)
plt.savefig('PairNeutrino.png')
plt.clf()

#plot#
plt.plot(plas[:,0], plas[:,1], label='1 MeV')
plt.plot(plas[:,0], plas[:,2], label='2 MeV')
plt.plot(plas[:,0], plas[:,3], label='3 MeV')
plt.plot(plas[:,0], plas[:,4], label='4 MeV')
plt.plot(plas[:,0], plas[:,5], label='5 MeV')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Production Rate (CODE)')
plt.title('Plasma Neutrino Spectra')
plt.grid(True)
plt.savefig('PlasmaNeutrino.png')
plt.clf()

#plot#
plt.plot(thermo[:,0], thermo[:,1]+thermo[:,2]+thermo[:,3]+thermo[:,4]+thermo[:,5])
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy Lost Rate (CODE)')
plt.title('Neutrino Energy Lost')
plt.grid(True)
plt.savefig('ThermoNeutrino.png')
plt.clf()
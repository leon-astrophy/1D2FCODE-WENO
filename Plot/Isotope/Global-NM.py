#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
mhe4=np.loadtxt('../../Outfile/Isotope/Star_WENO_MassHe4_1.dat',comments='"')
mc12=np.loadtxt('../../Outfile/Isotope/Star_WENO_MassC12_1.dat',comments='"')
mo16=np.loadtxt('../../Outfile/Isotope/Star_WENO_MassO16_1.dat',comments='"')
mne20=np.loadtxt('../../Outfile/Isotope/Star_WENO_MassNe20_1.dat',comments='"')
mmg24=np.loadtxt('../../Outfile/Isotope/Star_WENO_MassMg24_1.dat',comments='"')
msi28=np.loadtxt('../../Outfile/Isotope/Star_WENO_MassSi28_1.dat',comments='"')
mni56=np.loadtxt('../../Outfile/Isotope/Star_WENO_MassNi56_1.dat',comments='"')
ecap=np.loadtxt('../../Outfile/Isotope/Star_WENO_EcapNeutrinoEnergyLoss_1.dat',comments='"')

#plot#
plt.plot(mhe4[:,0], mhe4[:,1], label='He4')
plt.plot(mc12[:,0], mc12[:,1], label='C12')
plt.plot(mo16[:,0], mo16[:,1], label='O16')
plt.plot(mne20[:,0], mne20[:,1], label='Ne20')
plt.plot(mmg24[:,0], mmg24[:,1], label='Mg24')
plt.plot(msi28[:,0], msi28[:,1], label='Si28')
plt.plot(mni56[:,0], mni56[:,1], label='Ni56')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Mass (CODE)')
plt.title('Isotope Mass Versus Time')
plt.grid(True)
plt.savefig('isotope.png')
plt.clf()

#plot#
plt.plot(ecap[:,0], ecap[:,1])
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Energy Lost Rate (CODE)')
plt.title('Ecap-Neutrino Energy Lost Versus Time')
plt.grid(True)
plt.savefig('ecap.png')
plt.clf()
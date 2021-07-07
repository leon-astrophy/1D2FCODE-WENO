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
xc12=np.loadtxt('../../Outfile/Isotope/Star_WENO_XC12_1.dat',comments='"')
xo16=np.loadtxt('../../Outfile/Isotope/Star_WENO_XO16_1.dat',comments='"')
xhe4=np.loadtxt('../../Outfile/Isotope/Star_WENO_XHe4_1.dat',comments='"')
xmg24=np.loadtxt('../../Outfile/Isotope/Star_WENO_XMg24_1.dat',comments='"')
xsi28=np.loadtxt('../../Outfile/Isotope/Star_WENO_XSi28_1.dat',comments='"')
xne20=np.loadtxt('../../Outfile/Isotope/Star_WENO_XNe20_1.dat',comments='"')
xni56=np.loadtxt('../../Outfile/Isotope/Star_WENO_XNi56_1.dat',comments='"')

# get grid #
nmcount = int(np.size(xc12[:,0])/grid2)

#array#
nmtime = np.zeros(nmcount)

# load time array #
for i in range(0, nmcount):
    nmtime[i] = np.genfromtxt('../../Outfile/Isotope/Star_WENO_XC12_1.dat', skip_header = i * (grid2 + 2), max_rows = 1, usecols = (2))

#plot#
for i in range(0, nmcount):
   plt.plot(xhe4[i*grid2 + 5:(i+1)*grid2 - 5,0], xhe4[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('He$_{4}$')
plt.title('He$_{4}$ Profile')
plt.grid(True)
plt.savefig('he4.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(xc12[i*grid2 + 5:(i+1)*grid2 - 5,0], xc12[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('C$_{12}$')
plt.title('C$_{12}$ Profile')
plt.grid(True)
plt.savefig('c12.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(xo16[i*grid2 + 5:(i+1)*grid2 - 5,0], xo16[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('O$_{16}$')
plt.title('O$_{16}$ Profile')
plt.grid(True)
plt.savefig('o16.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(xne20[i*grid2 + 5:(i+1)*grid2 - 5,0], xne20[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Ne$_{20}$')
plt.title('Ne$_{20}$ Profile')
plt.grid(True)
plt.savefig('ne20.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(xmg24[i*grid2 + 5:(i+1)*grid2 - 5,0], xmg24[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Mg$_{24}$')
plt.title('Mg$_{24}$ Profile')
plt.grid(True)
plt.savefig('mg24.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(xsi28[i*grid2 + 5:(i+1)*grid2 - 5,0], xsi28[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Si$_{28}$')
plt.title('Si$_{28}$ Profile')
plt.grid(True)
plt.savefig('si28.png')
plt.clf()

#plot#
for i in range(0, nmcount):
   plt.plot(xni56[i*grid2 + 5:(i+1)*grid2 - 5,0], xni56[i*grid2 + 5:(i+1)*grid2 - 5,1], label = '%.3f' % nmtime[i])
plt.legend(loc="upper right")
plt.xlabel('Distance (CODE)')
plt.ylabel('Ni$_{56}$')
plt.title('Ni$_{56}$ Profile')
plt.grid(True)
plt.savefig('ni56.png')
plt.clf()

#prevent run-time error#
plt.close('all')
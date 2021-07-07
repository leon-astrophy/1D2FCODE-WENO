#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math

#loaddata#
lum=np.loadtxt('../../Outfile/Flame/Star_WENO_Luminosity_1.dat',comments='"')

#plot#
plt.plot(lum[:,0], lum[:,2], label='Flame')
plt.plot(lum[:,0], lum[:,3], label='Deton')
plt.plot(lum[:,0], lum[:,4], label='Burn')
plt.legend(loc="upper right")
plt.xlabel('Time (CODE)')
plt.ylabel('Luminosity (CODE)')
plt.title('Luminosity Versus Time')
plt.grid(True)
plt.savefig('lum.png')
plt.clf()

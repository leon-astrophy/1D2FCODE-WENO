#import required package#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

#loadeddata#
parameter=np.genfromtxt('../../Outfile/Parameter.dat')

#load parameter#
fileno = int(parameter[2,0])

for i in range (0,fileno+1):

   essential = np.loadtxt('../../Outfile/Tracer/Star_WENO_PPT_' + str(i) + '.dat', max_rows=1)
   time = essential[1]
   #load data#
   ppt = np.genfromtxt('../../Outfile/Tracer/Star_WENO_PPT_' + str(i) + '.dat', skip_header=2)

   # Generate Contour Plot #
   plt.plot(ppt[:,3], ppt[:,1], linestyle="--",label="Density")
   plt.xlabel('Distance')
   plt.ylabel('Density')
   plt.grid('True')
   plt.title('Tracer Density Profile For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.savefig('PPT-Tracer-Density'+str(i))
   plt.clf()

   # Generate Contour Plot #
   plt.plot(ppt[:, 3], ppt[:, 2], linestyle="--", label="Temperature")
   plt.xlabel('Distance')
   plt.ylabel('Temperature')
   plt.grid('True')
   plt.title('Tracer Temperature Profile For Time = ' + str("{:.3f}".format(time)), y=1.08, fontsize=11)
   plt.tight_layout()
   plt.savefig('PPT-Tracer-Temperature'+str(i))
   plt.clf()

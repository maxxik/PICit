import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import ticker
import sys
from cycler import cycler

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


matplotlib.rc("font", family="serif")
matplotlib.rc("font", size=14)

plt.rc('axes', prop_cycle=(cycler(linestyle=['-', '--', '-.','dotted',])*
cycler('color', ['r', 'navy','g','orange','y','c','m','k','gray','b','olive','brown'])))

name=sys.argv[1]
csnum=int(sys.argv[2])
Ntarget=int(sys.argv[3])

print('PICit test_cross_section() figure generation...')

A=np.loadtxt(name)
print(A.shape)

en=A[::10,0]




print ('Plotting cross sections...\n ')
plt.figure(figsize=cm2inch(15,10))
plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0, wspace=0.1, hspace=0.1)
ax1 = plt.subplot2grid((1,1), (0,0))
for i in range(csnum):
	ax1.plot(en,A[::10,Ntarget+i+1],lw=0.7,label=r"Process {}".format(i+1))



ax1.grid(True)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel(r'$\epsilon$ [eV]')
ax1.set_ylabel(r'$\sigma$ [m$^2$]')
ax1.set_xlim(en[1],en[-1])
ax1.set_ylim(1e-26,)

ax1.legend(loc=1, ncol=2, bbox_to_anchor=(0.5, -0.15, 1.2, 1.2))

plt.savefig(name+'.png',format="png",bbox_inches='tight', dpi=600)

plt.close()

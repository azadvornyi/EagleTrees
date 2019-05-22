
import pylab as pl
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))


keep_these = np.where(t_q != 14.134423953091941) #
and_these = np.where((t_i<=t_q))
q_before_i = np.where((t_i>=t_q))
unquenched = np.where(t_q == 14.134423953091941)

a =set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
c = np.array(list(c))
font = {'family' : 'monospace',
        'size'   : '12'}
plt.rc('font', **font)
#pl.scatter(t_i,t_q,lw=0.5,c='k',edgecolor='w')  #overlaying the sample points
fig = plt.figure(figsize=[6,6])

ax = fig.add_subplot(111)
ax.set_aspect('equal')                           #sub-plot area 2 out of 2

heatmap = ax.hexbin(t_i[c],t_q[c],C=None,gridsize=15,bins=None,mincnt=0,cmap='Blues')
ax.plot(np.linspace(0,14,10),np.linspace(0,14,10), c='red')
#cbar = fig.colorbar(heatmap, ax=ax)

cbar = plt.colorbar(heatmap,fraction=0.046, pad=0.04)
cbar.set_label(r'$\mathtt{N}$', rotation=90)
#ax2 = ax.twiny()
ax.set_aspect('equal')
minorLocatory = AutoMinorLocator()
minorLocatorx = AutoMinorLocator()


#cbar.minorticks_on()

ax.text(8.4,9.2, r'$ \mathtt{t_{quench} = t_{infall}}  $',
         rotation=45,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax.xaxis.set_minor_locator(minorLocatorx)
ax.yaxis.set_minor_locator(minorLocatory)
ax.tick_params(axis='both', direction='in', length=6, width=1.1,which = 'major', colors='black',
               grid_color='black', grid_alpha=0.4)
ax.tick_params(which='minor', length=4, color='black',direction='in')
ax.axis('image')                            #necessary for correct aspect ratio
ax.set_xlim(2,14)
ax.set_ylim(2,14)
ax.set_ylabel(r'$\mathtt{t_{quench}[Gyr]}$')
ax.set_xlabel(r'$\mathtt{t_{infall}[Gyr]}$')
plt.savefig('hexbins.pdf', format ='pdf')
plt.show()

import pylab as pl
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from pylab import *
M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))


keep_these = np.where(t_q != 14.134423953091941) #
and_these = np.where((t_i<=t_q))
q_before_i = np.where((t_i>=t_q))
unquenched = np.where(t_q == 14.134423953091941)
nonzero = np.where(sfr!=0)[0]
a =set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
c = np.array(list(c))
font = {'family' : 'monospace',
        'size'   : '12'}
#plt.rc('text', usetex=False)
#pl.scatter(t_i,t_q,lw=0.5,c='k',edgecolor='w')  #overlaying the sample points
fig = plt.figure(figsize=[6.5,6])

ax = fig.add_subplot(111)
#ax.set_aspect('equal', 'box')                           #sub-plot area 2 out of 2

heatmap = ax.hexbin(np.log10(m_s[nonzero]),np.log10(sfr[nonzero]),C=None,gridsize=17,bins=None,mincnt=0,cmap='Blues')
#ax.scatter(np.log10(m_s),(sfr)*1e9)
ax.hlines(-11,8.9,max(np.log10(m_s)), colors="red")
ax.set_xlim(8.9,11.5)
#ax.plot(np.linspace(0,14,10),np.linspace(0,14,10), c='red')
#cbar = fig.colorbar(heatmap, ax=ax)
#ax.set_ylim(0, 10**(-11))
cbar = plt.colorbar(heatmap,fraction=0.046, pad=0.04)
cbar.set_label(r'$\rm N$', rotation=90)
#ax.set_aspect('equal', 'box')
#ax.set_aspect('equal')
# minorLocatory = AutoMinorLocator()
# minorLocatorx = AutoMinorLocator()
# minorLocatory2 = AutoMinorLocator()
# minorLocatorx2 = AutoMinorLocator()
# minorLocatory3 = AutoMinorLocator()
# minorLocatorx3 = AutoMinorLocator()
# #cbar.minorticks_on()
#
# # ax.text(9,5,-11.8, r'$ \mathtt{SSFRcut}$',
# #          rotation=0,
# #          horizontalalignment='center',
# #          verticalalignment='top',
# #          multialignment='center')
#
# ax.xaxis.set_minor_locator(minorLocatorx)
# ax.yaxis.set_minor_locator(minorLocatory)
#
# ax.tick_params(axis='both', direction='in', length=6, width=1.1,which = 'major', colors='black',
#                grid_color='black', grid_alpha=0.4)
# ax.tick_params(which='minor', length=4, color='black',direction='in')

#ax.axis('image')                            #necessary for correct aspect ratio
ax.set_ylabel(r'$\rm log_{10}(SSFR/yr^{-1})$')
ax.set_xlabel(r'$\rm log_{10}( \it M_{\star} \rm /M_{\odot})$')
plt.savefig('plots_to_combine/red_valley.pdf', format ='pdf', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.show()
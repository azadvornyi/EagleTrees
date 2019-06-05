import matplotlib as mpl

mpl.use('TkAgg')
import matplotlib.pyplot as plt
import time
from scipy import interpolate
from sys import argv
import os
import pickle
import numpy as np
from scipy import stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys

Rperi_gid, Rperi_satn, Rperi,M_h, m_s, t_i, t_q, sfr, m_sat_max, v_peri, mhost = np.loadtxt("Rperi", unpack=True, usecols=(0,1,2, 3, 4, 5, 6,7, 8,9,10))

#p_ram = (v_peri*(mhost**(2/3))) / (Rperi**2)
p_ram = v_peri
# maxRperi = np.where(Rperi>0.9)[0]
# print(Rperi_gid[maxRperi], Rperi_satn[maxRperi], Rperi[maxRperi])
# print(max(Rperi))

# plt.hist(np.log10(p_ram), bins=30)
#
# plt.ylabel('N')
# plt.xlabel(r'$p{ram} \ [(km/s)^2 M_{\star}^{2/3} cMpc^{-2}]$')
# plt.savefig('p_ram_distribution.pdf',format='pdf')
# plt.show()
# sys.exit()


logmassmin = 9
logmassmax = 12
dex = 0.25
bin_num = (logmassmax - logmassmin) / dex
bin_range = np.linspace(9, 12, bin_num + 1)

percentile25 = np.percentile(p_ram, 25)
percentile50 = np.percentile(p_ram, 50)
percentile75 = np.percentile(p_ram, 75)


#
# plt.hist(np.sqrt(p_ram), bins=30, color = 'grey')
#
# plt.ylabel('N')
# plt.xlabel(r'$v_{peri} \ [km/s]$')
#
# plt.vlines(np.sqrt(percentile25), 0, 125, color='black')
#
# plt.vlines(np.sqrt(percentile50), 0, 125, color='black')
#
# plt.vlines(np.sqrt(percentile75), 0, 125, color='black')
# plt.title("v_peri  distribution")
# #plt.legend()
# #plt.savefig('Rperi_distribution.pdf',format='pdf')
# #plt.savefig('plots_to_combine/v2_distribution.jpg',format='jpg',dpi=1200, pad_inches=0.1, bbox_inches='tight')
# plt.show()
# sys.exit()



#
# small_split = np.where(p_ram < percentile25)[0]
# large_split = np.where(p_ram > percentile75)[0]

small_split = np.where(p_ram <= percentile50)[0]
large_split = np.where(p_ram>percentile50)[0]

print(percentile50)

quench_timescale_s = t_q[small_split] - t_i[small_split]
quench_timescale_l = t_q[large_split] - t_i[large_split]
# bin_means, bin_edges, binnumber = stats.binned_statistic( np.log10(m_s[c]),quench_timescale, statistic='median', bins=[np.log10(1e9),
#                             np.log10(10**9.5) ,np.log10(1e10),np.log10(10**10.5) ,np.log10(1e11), np.log10(10**11.5) ,np.log10(1e12)])

bin_means_s, bin_edges_s, binnumber_s = stats.binned_statistic(np.log10(m_s[small_split]), quench_timescale_s,
                                                               statistic='median', bins=bin_range)
bin_means_l, bin_edges_l, binnumber_l = stats.binned_statistic(np.log10(m_s[large_split]), quench_timescale_l,
                                                               statistic='median', bins=bin_range)
bin_width_s = (bin_edges_s[1] - bin_edges_s[0])
bin_centers_s = bin_edges_s[1:] - bin_width_s / 2

bin_width_l = (bin_edges_l[1] - bin_edges_l[0])
bin_centers_l = bin_edges_l[1:] - bin_width_l / 2

# fig = plt.figure(figsize=(8, 6))
# ax = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan=2)
#
# ax.hlines(bin_means_l, bin_edges_l[:-1] + 0.01, bin_edges_l[1:] + 0.01, colors='#66a5ad', lw=1,label="upper")
# ax.hlines(bin_means_s - 0.01, bin_edges_s[:-1] - 0.01, bin_edges_s[1:] - 0.01, colors='#07575b', lw=1, label="lower")
# ax.legend()
# #plt.hist(Rperi, bins = 28)
# plt.show()

bin_means75_s, bin_edges75_s, binnumber75_s = stats.binned_statistic(np.log10(m_s[small_split]), quench_timescale_s,
                                                                     statistic=lambda x75_s: np.percentile(x75_s, 75),
                                                                     bins=bin_range)
bin_means25_s, bin_edges25_s, binnumber25_s = stats.binned_statistic(np.log10(m_s[small_split]), quench_timescale_s,
                                                                     statistic=lambda x25_s: np.percentile(x25_s, 25),
                                                                     bins=bin_range)

bin_means75_l, bin_edges75_l, binnumber75_l = stats.binned_statistic(np.log10(m_s[large_split]), quench_timescale_l,
                                                                     statistic=lambda x75_l: np.percentile(x75_l, 75),
                                                                     bins=bin_range)
bin_means25_l, bin_edges25_l, binnumber25_l = stats.binned_statistic(np.log10(m_s[large_split]), quench_timescale_l,
                                                                     statistic=lambda x25_l: np.percentile(x25_l, 25),
                                                                     bins=bin_range)
# print(type(bin_means))

# x = np.log10(m_s[c])
# bins = np.array([np.log10(1e9),np.log(10**9.5) ,np.log(1e10),np.log(10**10.5) ,np.log(1e11), np.log(10**11.5) ,np.log(1e12)])
# digitized = np.digitize(x, bins)
# print(digitized)
# print(x)
# print(np.log10(10))
# bin_means = [x[digitized == i].mean() for i in range(1, len(bins))]
# print(bin_means)
font = {'family': 'monospace',
        'size': '12'}
header = {'family': 'monospace',
          'size': '15'}
small_font = {'family': 'monospace',
              'size': '10'}
plt.rc('font', **font)



# plot it
fig = plt.figure(figsize=(8, 6))
ax = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan=2)

ax5 = plt.subplot2grid((3, 2), (2, 0), colspan=2, rowspan=1)
ax5.hist(np.log10(m_s[small_split]), bins=bin_edges_s, align='mid', color='#07575b',
         label=r'$\mathtt{Q_1+Q_2}$', edgecolor='#003b46', linewidth=1.2,alpha = 0.65)
ax5.hist(np.log10(m_s[large_split]), bins=bin_edges_s, align='mid', label=r'$\mathtt{Q_3+Q_4}$',
         color='#66a5ad', edgecolor='#003b46', linewidth=1.2, alpha = 0.65)

# plt.yscale('log')
# plt.xlabel("M_sat [M_sun]")
# plt.ylabel("N")
minorLocatorx = MultipleLocator(0.25)
minorLocatory = AutoMinorLocator()

ax.xaxis.set_minor_locator(minorLocatorx)
ax.yaxis.set_minor_locator(minorLocatory)
ax.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
               grid_color='black', grid_alpha=0.4, labelbottom=False)
ax.tick_params(which='minor', length=4, color='black', direction='in')
ax.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
               grid_color='black', grid_alpha=0.4)
ax2 = ax.twiny()

ax2.xaxis.set_minor_locator(minorLocatorx)
ax2.yaxis.set_minor_locator(minorLocatory)
ax2.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4, labeltop=False, labelright=False)
ax2.tick_params(which='minor', length=4, color='black', direction='in')
ax2.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax2.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
                grid_color='black', grid_alpha=0.4)

ax3 = ax.twinx()
minorLocatory = AutoMinorLocator()
minorLocatorx = MultipleLocator(0.25)
ax3.xaxis.set_minor_locator(minorLocatorx)
ax3.yaxis.set_minor_locator(minorLocatory)
ax3.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4)
ax3.tick_params(which='minor', length=4, color='black', direction='in')
ax3.tick_params(axis='y', which='minor', length=4, color='black', direction='in', labelright=False)
ax3.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
                grid_color='black', grid_alpha=0.4)
ax3.set_ylim(0, 8)
ax2.set_xlim(8.5, 12)
ax3.set_yticklabels((' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '))

print(bin_means75_s)
print(bin_means25_s)
print(bin_means_s)

ax.hlines(bin_means_s - 0.01, bin_edges_s[:-1] - 0.01, bin_edges_s[1:] - 0.01, colors='#07575b', lw=1)

ax.vlines(bin_centers_s - 0.01, bin_means25_s, bin_means75_s, colors='#07575b', lw=1, )
# ax.hlines(bin_means_s[:-3 or None]-0.01, bin_edges_s[:-1]-0.01, bin_edges_s[1:]-0.01, colors='g', lw=1.5)
# ax.vlines(bin_centers_s[:-3 or None]-0.01,bin_means25_s[:-3 or None],bin_means75_s[:-3 or None],colors='g', lw=1.5,)

ax.hlines(bin_means_l, bin_edges_l[:-1] + 0.01, bin_edges_l[1:] + 0.01, colors='#66a5ad', lw=1)

ax.vlines(bin_centers_l + 0.01, bin_means25_l, bin_means75_l, colors='#66a5ad', lw=1, )

ax.errorbar(1, 1, 1, 1, label=r'$\mathtt{Q_1+Q_2}$', c='#07575b', lw=1)

ax.errorbar(1, 1, 1, 1, label=r'$\mathtt{Q_3+Q_4}$', c='#66a5ad', lw=1)

ax.set_ylabel(r"$\mathtt{ \tau_q [Gyr]}$")
# plt.plot((binnumber - 0.5) * bin_width, x_pdf, 'g.', alpha=0.5)
ax5.set_xlabel(r"$\mathtt{log_{10}(M_{\bigstar}/M_{\odot}}) $")
# plt.xscale('log')
ax.set_ylim(0, 8)
ax2.set_xlim(8.9, 11.5)
ax.set_xlim(8.9, 11.5)
ax5.set_xlim(8.9, 11.5)
ax5.set_ylabel(r'$\mathtt{N}$')
minorLocatory5 = AutoMinorLocator()
minorLocatorx5 = MultipleLocator(0.25)
ax5.xaxis.set_minor_locator(minorLocatorx5)
ax5.yaxis.set_minor_locator(minorLocatory5)
ax5.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4)
ax5.tick_params(which='minor', length=4, color='black', direction='in')
ax5.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax5.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
                grid_color='black', grid_alpha=0.4)
ax.legend(frameon=False)
ax.set_title(r'$\mathtt{QTD \ \ split \ \  by \ \ v^2_{peri} }$')
#ax.set_title(r'$\mathtt{QTD \ \ split \ \  by \ \ P_{ram} }$')
ax5.legend(frameon=False)
plt.subplots_adjust(hspace=0, wspace=0.33)

#plt.savefig('qtd_1234_v2_peri.pdf', format='pdf', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.savefig('plots_to_combine/qtd_1234_v2_peri.jpg', format='jpg', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.show()

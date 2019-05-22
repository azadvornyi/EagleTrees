import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504.txt", unpack=True, usecols=(2, 3, 4, 5, 6))
# M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RecalL0100N0752.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))

n15 = np.where(M_h > 10 ** 15)
print(len(M_h[n15]))
print(max(M_h))
logmassmin = 9
logmassmax = 12
dex = 0.25
bin_num = (logmassmax - logmassmin) / dex
bin_range = np.linspace(9, 12, bin_num + 1)

keep_these = np.where(t_q != 14.134423953091941)  #
and_these = np.where((t_i <= t_q))
q_before_i = np.where((t_i >= t_q))
unquenched = np.where(t_q == 14.134423953091941)
# print(type(and_these))
# print(len(keep_these))
# print(and_these[0])
nuq = len(m_s) - len(m_s[keep_these]) + len(m_s[and_these])
print(keep_these)
# print(m_s[keep_these])
a = set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
c = np.array(list(c))

ratio = m_s[c]/M_h[c]

# plt.hist(np.log10(ratio), bins=30)
#
# plt.ylabel('N')
# plt.xlabel(r'$log_{10}(M_{*}/M_{host})$')
# plt.savefig('log10_ratio_distribution.pdf',format='pdf')
# plt.show()
# sys.exit()
# print("quenched after infall: {0} \nquenched before infall: {1}\n unquenched: {2}\n total {3}".format(len(c),
#                                                                                 len(m_s[q_before_i]),len(m_s[unquenched]),len(m_s)))
percentile25 = np.percentile(ratio, 25)
percentile50 = np.percentile(ratio, 50)
percentile75 = np.percentile(ratio, 75)

small_split = np.where(ratio < percentile25)[0]
# small_split = np.where(ratio <= percentile50)[0]
# large_split = np.where(ratio>percentile50)[0]
large_split = np.where(ratio > percentile75)[0]

quench_timescale = t_q[c] - t_i[c]
quench_timescale_s = t_q[small_split] - t_i[small_split]
quench_timescale_l = t_q[large_split] - t_i[large_split]
# bin_means, bin_edges, binnumber = stats.binned_statistic( np.log10(m_s[c]),quench_timescale, statistic='median', bins=[np.log10(1e9),
#                             np.log10(10**9.5) ,np.log10(1e10),np.log10(10**10.5) ,np.log10(1e11), np.log10(10**11.5) ,np.log10(1e12)])

bin_means_s, bin_edges_s, binnumber_s = stats.binned_statistic(np.log10(m_s[small_split]), quench_timescale_s,
                                                               statistic='median', bins=bin_range)
bin_means_l, bin_edges_l, binnumber_l = stats.binned_statistic(np.log10(m_s[large_split]), quench_timescale_l,
                                                               statistic='median', bins=bin_range)
# print(bin_means,bin_edges, binnumber)
bin_width_s = (bin_edges_s[1] - bin_edges_s[0])
bin_centers_s = bin_edges_s[1:] - bin_width_s / 2

bin_width_l = (bin_edges_l[1] - bin_edges_l[0])
bin_centers_l = bin_edges_l[1:] - bin_width_l / 2

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
         label=r'$\mathtt{Q_1}$', edgecolor='#003b46', linewidth=1.2, alpha = 0.65)
ax5.hist(np.log10(m_s[large_split]), bins=bin_edges_s, align='mid', label=r'$\mathtt{Q_4}$',
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

ax.errorbar(1, 1, 1, 1, label=r'$\mathtt{Q_1}$', c='#07575b', lw=1)

ax.errorbar(1, 1, 1, 1, label=r'$\mathtt{Q_4}$', c='#66a5ad', lw=1)

ax.set_ylabel(r"$\mathtt{ \tau_q [Gyr]}$")
# plt.plot((binnumber - 0.5) * bin_width, x_pdf, 'g.', alpha=0.5)
ax5.set_xlabel(r"$\mathtt{log_{10}(M_{\bigstar}/M_{\odot}}) $")
# plt.xscale('log')
ax.set_ylim(0, 8)
ax.set_xlim(8.9, 12)
ax2.set_xlim(8.9, 12)
ax5.set_xlim(8.9, 12)
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
ax5.legend(frameon=False)
ax.set_label(r'$\mathtt{Split by M_{*}/M_{host}}$')
plt.subplots_adjust(hspace=0, wspace=0.33)
#plt.savefig('qtd_df_Q1_Q4 .pdf', format='pdf', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.font_manager
M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))
#M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RecalL0025N0752.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))


logmassmin = 9
logmassmax = 12
dex = 0.25
bin_num = (logmassmax - logmassmin)/dex + 1
bin_range = np.linspace(9, 12, bin_num)

keep_these = np.where(t_q != 14.134423953091941) #
and_these = np.where((t_i<=t_q))
q_before_i = np.where((t_i>=t_q))
unquenched = np.where(t_q == 14.134423953091941)
print(type(and_these))
# print(len(keep_these))
#print(and_these[0])
nuq = len(m_s)- len(m_s[keep_these]) + len(m_s[and_these])

#print(m_s[keep_these])
a =set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
# c = np.array(list(c))


m9 = set((np.where((m_s >= 10**9) & (m_s <10**9.5))[0]))
m95 = set((np.where((m_s >= 10**9.5) & (m_s <10**10))[0]))
m10 = set((np.where((m_s >= 10**10) & (m_s <10**10.5))[0]))
m105 = set((np.where((m_s >= 10**10.5) & (m_s <10**11))[0]))
m11 = set((np.where((m_s >= 10**11))[0]))


c9 = c.intersection(m9)
c95 = c.intersection(m95)
c10= c.intersection(m10)
c105= c.intersection(m105)
c11= c.intersection(m11)

c9 = np.array(list(c9))
c95 = np.array(list(c95))
c10 = np.array(list(c10))
c105 = np.array(list(c105))
c11 = np.array(list(c11))

c = np.array(list(c))
import matplotlib
matplotlib.matplotlib_fname()

#font = {'family' : 'monospace',
#        'size'   : '12'}

font = {'size'   : '12'}
plt.rc('text', usetex=False)
print("quenched after infall: {0} \nquenched before infall: {1}\n unquenched: {2}\n total {3}".format(len(c),
                                                                                len(m_s[q_before_i]),len(m_s[unquenched]),len(m_s)))

quench_timescale = t_q[c] - t_i[c]
# bin_means, bin_edges, binnumber = stats.binned_statistic( np.log10(m_s[c]),quench_timescale, statistic='median', bins=[np.log10(1e9),
#                             np.log10(10**9.5) ,np.log10(1e10),np.log10(10**10.5) ,np.log10(1e11), np.log10(10**11.5) ,np.log10(1e12)])

bin_means, bin_edges, binnumber = stats.binned_statistic( np.log10(m_s[c]),quench_timescale, statistic='median', bins=bin_range)

print(bin_means,bin_edges, binnumber)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2

#plt.rc('text', usetex=False)

# x = np.log10(m_s[c])
# bins = np.array([np.log10(1e9),np.log(10**9.5) ,np.log(1e10),np.log(10**10.5) ,np.log(1e11), np.log(10**11.5) ,np.log(1e12)])
# digitized = np.digitize(x, bins)
# print(digitized)
# print(x)
# print(np.log10(10))
# bin_means = [x[digitized == i].mean() for i in range(1, len(bins))]
# print(bin_means)
fig = plt.figure(figsize=[8,8])
plt.subplot(321)
plt.scatter(t_q[c]-t_i[c],t_i[c], s= 2, c= "blue", label = r"$\mathrm{9 < log_{10}(M_{\star}/M_{\odot})<11}$")
plt.ylabel(r"$\tau_q [Gyr]$")
plt.xlabel(r"$t_i \ [Gyr]$")
plt.ylim(-0.5, 14.2)
plt.xlim(-0.5, 14.2)

plt.legend()

plt.subplot(322)
plt.scatter(t_q[c9]-t_i[c9],t_i[c9], s= 2, c= "blue", label = r"$ \rm 9 < \rm log_{10}( \it M_{\star} \rm /M_{\odot})<9.5$")
plt.ylabel(r"$\tau_q [Gyr]$")
plt.xlabel(r"$t_i \ [Gyr]$")
plt.ylim(-0.5, 14.2)
plt.xlim(-0.5, 14.2)

plt.legend()
#================
plt.subplot(323)
plt.scatter(t_q[c95]-t_i[c95],t_i[c95], s= 2, c= "blue", label = r"$9.5 < log_{10}(M_{\star}/M_{\odot})<10$")
plt.ylabel(r"$\tau_q [Gyr]$")
plt.xlabel(r"$t_i \ [Gyr]$")
plt.ylim(-0.5, 14.2)
plt.xlim(-0.5, 14.2)
plt.legend()

plt.subplot(324)
plt.scatter(t_q[c10]-t_i[c10],t_i[c10], s= 2, c= "blue", label = r"$10 < log_{10}(M_{\star}/M_{\odot})<10.5$")
plt.ylabel(r"$\tau_q [Gyr]$")
plt.xlabel(r"$t_i \ [Gyr]$")
plt.ylim(-0.5, 14.2)
plt.xlim(-0.5, 14.2)
plt.legend()
#==============plot
plt.subplot(325)
plt.scatter(t_q[c105]-t_i[c105],t_i[c105], s= 2, c= "blue", label = r"$10.5 < log_{10}(M_{\star}/M_{\odot})<11$")
plt.ylabel(r"$\tau_q [Gyr]$")
plt.xlabel(r"$t_i \ [Gyr]$")
plt.ylim(-0.5, 14.2)
plt.xlim(-0.5, 14.2)
plt.legend()


plt.subplot(326)
plt.scatter(t_q[c11]-t_i[c11],t_i[c11], s= 2, c= "blue", label = r"$11 < log_{10}(M_{\star}/M_{\odot})$")
plt.ylabel(r"$\tau_q [Gyr]$")
plt.xlabel(r"$t_i \ [Gyr]$")
plt.ylim(-0.5, 14.2)
plt.xlim(-0.5, 14.2)
plt.legend()
plt.tight_layout()
#plt.savefig('3plots_100.pdf', format ='pdf')
#plt.savefig('plots_to_combine/tauq_ti_dex_split.jpg', format='jpg', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504_old.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))




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
c = np.array(list(c))

print("quenched after infall: {0} \nquenched before infall: {1}\n unquenched: {2}\n total {3}".format(len(c),
                                                                                len(m_s[q_before_i]),len(m_s[unquenched]),len(m_s)))

quench_timescale = t_q[c] - t_i[c]
bin_means, bin_edges, binnumber = stats.binned_statistic( np.log10(m_s[c]),quench_timescale, statistic='median', bins=[np.log10(1e9),
                            np.log10(10**9.5) ,np.log10(1e10),np.log10(10**10.5) ,np.log10(1e11), np.log10(10**11.5) ,np.log10(1e12)])

print(bin_means,bin_edges, binnumber)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2



# x = np.log10(m_s[c])
# bins = np.array([np.log10(1e9),np.log(10**9.5) ,np.log(1e10),np.log(10**10.5) ,np.log(1e11), np.log(10**11.5) ,np.log(1e12)])
# digitized = np.digitize(x, bins)
# print(digitized)
# print(x)
# print(np.log10(10))
# bin_means = [x[digitized == i].mean() for i in range(1, len(bins))]
# print(bin_means)
fig = plt.figure(figsize=[9,8])
plt.subplot(321)
plt.scatter(t_i[c],t_q[c], s= 2, c= "blue")
plt.ylabel("t_q [Gyr]")
plt.xlabel("t_i [Gyr]")

plt.subplot(323)
plt.hist(t_i[c],bins=28)
plt.xlabel("t_i [Gyr]")
plt.ylabel("N")
#================
plt.subplot(322)
plt.scatter(m_s[c], t_q[c]-t_i[c],s= 2, c= "blue")
plt.xscale('log')
plt.ylabel("t_q - t_i [Gyr]")
plt.xlabel("M_sat [M_sun]")

plt.subplot(324)
plt.hist(np.log(m_s[c]), bins=28)
plt.xscale('log')
plt.xlabel("M_sat [M_sun]")
plt.ylabel("N")
#==============
plt.subplot(325)
plt.scatter(t_i[c], sfr[c]*1e9, s= 2, c= "blue")
plt.ylabel("ssfr(z=0)  Gyr^-1 ")
plt.xlabel("t_i [Gyr]")
# plt.yscale('log')
# plt.ylim(1e-9,1)
plt.tight_layout()

plt.subplot(326)
plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=1,
           label='binned statistic of data')
plt.ylabel("Quenching timescale [Gyr]")
#plt.plot((binnumber - 0.5) * bin_width, x_pdf, 'g.', alpha=0.5)
plt.xlabel("M_sat [M_sun]")
#plt.xscale('log')

plt.savefig('3plots.pdf', format ='pdf')
plt.show()
import numpy as np
import matplotlib.pyplot as plt

M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_100.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))




keep_these = np.where(t_q != 14.134423953091941)
# print(len(keep_these))
nuq = len(m_s)- len(m_s[keep_these])
# print(m_s[keep_these])

print("quenched: {0} \ntotal: {1}".format(nuq, len(m_s)))
fig = plt.figure(figsize=[9,8])
plt.subplot(321)
plt.scatter(t_i[keep_these],t_q[keep_these], s= 2, c= "blue")
plt.ylabel("t_q [Gyr]")
plt.xlabel("t_i [Gyr]")

plt.subplot(323)
plt.hist(t_i[keep_these],bins=28)
plt.xlabel("t_i [Gyr]")
plt.ylabel("N")
#================
plt.subplot(322)
plt.scatter(m_s[keep_these], t_q[keep_these]-t_i[keep_these],s= 2, c= "blue")
plt.xscale('log')
plt.ylabel("t_q - t_i [Gyr]")
plt.xlabel("M_sat [M_sun]")

plt.subplot(324)
plt.hist(np.log(m_s[keep_these]), bins=28)
plt.xscale('log')
plt.xlabel("M_sat [M_sun]")
plt.ylabel("N")
#==============
plt.subplot(325)
plt.scatter(t_i[keep_these], sfr[keep_these]*1e9, s= 2, c= "blue")
plt.ylabel("ssfr(z=0)  Gyr^-1 ")
plt.xlabel("t_i [Gyr]")
# plt.yscale('log')
# plt.ylim(1e-9,1)
plt.tight_layout()
plt.savefig('3plots.pdf', format ='pdf')
plt.show()
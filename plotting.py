import numpy as np
import matplotlib.pyplot as plt

M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot.txt", unpack="True", usecols=(0, 1, 2, 3, 4,))

fig = plt.figure(figsize=[9,8])
plt.subplot(321)
plt.scatter(t_i,t_q, s= 2, c= "blue")
plt.ylabel("t_q [Gyr]")
plt.xlabel("t_i [Gyr]")

plt.subplot(323)
plt.hist(t_i,bins=28)
plt.xlabel("t_i [Gyr]")
plt.ylabel("N")
#================
plt.subplot(322)
plt.scatter(m_s, t_q-t_i,s= 2, c= "blue")
plt.xscale('log')
plt.ylabel("t_q - t_i [Gyr]")
plt.xlabel("M_sat [M_sun]")

plt.subplot(324)
plt.hist(m_s,bins=56)
plt.xscale('log')
plt.xlabel("M_sat [M_sun]")
plt.ylabel("N")
#==============
plt.subplot(325)
plt.scatter(t_i, sfr*1e9, s= 2, c= "blue")
plt.ylabel("sfr(z=0)  M_s/Gyr ")
plt.xlabel("t_i [Gyr]")
plt.tight_layout()
plt.savefig('3plots.pdf', format ='pdf')
plt.show()
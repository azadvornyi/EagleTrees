import matplotlib.pyplot as plt
import time
from scipy import interpolate
from sys import argv
import os
import scipy as sc
import pickle
import numpy as np
from scipy.signal import argrelextrema
import sys
from astropy.cosmology import FlatLambdaCDM
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.lines as mlines
import astropy.units as u
from astropy.cosmology import z_at_value
from collections import deque,Counter
from bisect import insort, bisect_left
from itertools import islice

host_gid, sat_num, M_h, m_s_max, t_i, t_q, tau_q, sfr = np.loadtxt("Data_plot_100_with_max_mass", unpack=True,
                                                                   usecols=(0, 1, 2, 3, 4, 5, 6, 7,))

unique_hosts, unique_hosts_index = np.unique(host_gid, return_index=True)



plt.rc('text', usetex=False)


dex = 0.25
dex_range = np.linspace(13,14.5, 7)

def smooth_peaks(dex_edge, dex):

    h325 = np.where((np.log10(M_h[unique_hosts_index]) > dex_edge) & (np.log10(M_h[unique_hosts_index]) < dex_edge + dex))[0]
    hosts325 = host_gid[unique_hosts_index][h325]
    #print(hosts325)
    sats325_index = np.array([])
    tau_q325 = []
    m_s_max325 = []
    for i in hosts325:
        holder = np.where(host_gid == i)[0]
        tau_q325 = np.append(tau_q325, tau_q[holder])
        m_s_max325 = np.append(m_s_max325, m_s_max[holder])


    tup325 = []
    for j in range(len(m_s_max325)):
        t = (m_s_max325[j],tau_q325[j])
        tup325.append(t)

    #print(tup325)
    #tup325 = zip(m_s_max325, tau_q325)

    # print(tup325)
    sorted_by_second = sorted(tup325, key=lambda tup: tup[0])
    # print(sorted_by_second)
    # print(sorted_by_second[2][1])

    W = 201
    m = int(W/2- 0.5)
    seq = iter(sorted_by_second)
    s = [item for item in islice(seq,W)]
    d = deque(s)
    median_tau = lambda: sorted(s,key=lambda tup: tup[1])[m][1]
    median_m = lambda : s[m][0]
    medians_tau = [median_tau()]
    medians_m = [median_m()]
    for item in seq:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        medians_tau.append(median_tau())
        medians_m.append(median_m())
    return medians_tau, medians_m

data_t = []
data_m = []
for i in dex_range:
    data_t.append(smooth_peaks(i,dex)[0])
    data_m.append(smooth_peaks(i, dex)[1])

print(data_t[0][0])
h = np.where((np.log10(M_h[unique_hosts_index]) > 13) & (np.log10(M_h[unique_hosts_index]) < 13 + dex))[0]
hosts = host_gid[unique_hosts_index][h]
#print(hosts325)
sats_index = np.array([])
tau_q1 = []
m_s_max1 = []
for i in hosts:
    holder = np.where(host_gid == i)[0]
    tau_q1 = np.append(tau_q1, tau_q[holder])
    m_s_max1 = np.append(m_s_max1, m_s_max[holder])


tup = []
for j in range(len(m_s_max1)):
    t = (m_s_max1[j],tau_q1[j])
    tup.append(t)

#print(tup325)
#tup325 = zip(m_s_max325, tau_q325)

#print(tup325)
sorted_by_second = sorted(tup, key=lambda tup: tup[0])[:25]
tau = sorted(sorted_by_second, key=lambda tup: tup[1])[11]
print(tau)
fig = plt.figure(figsize=[8,8])
alpha = 0.7
plt.plot(np.log10(data_m[0]),data_t[0], label = r'$\rm log_{10}(M_h/M_{\odot}) =(13.00 ,13.25) $',alpha=alpha,c='red')
plt.plot(np.log10(data_m[1]),data_t[1], label = r'$\rm log_{10}(M_h/M_{\odot}) =(13.25 ,13.50) $',alpha=alpha,c='orange')
plt.plot(np.log10(data_m[2]),data_t[2], label = r'$\rm log_{10}(M_h/M_{\odot}) =(13.50 ,13.75) $',alpha=alpha,c='#666666')
plt.plot(np.log10(data_m[3]),data_t[3], label = r'$\rm log_{10}(M_h/M_{\odot}) =(13.75 ,14.00) $',alpha=alpha,c='green')
plt.plot(np.log10(data_m[4]),data_t[4], label = r'$\rm log_{10}(M_h/M_{\odot}) =(14.00 ,14.25) $',alpha=alpha,c='blue')
plt.plot(np.log10(data_m[5]),data_t[5], label = r'$\rm log_{10}(M_h/M_{\odot}) =(14.25 ,14.50) $',alpha=alpha,c='#CD00CD')
plt.plot(np.log10(data_m[6]),data_t[6], label = r'$\rm log_{10}(M_h/M_{\odot}) =(14.50 ,14.75) $',alpha=alpha,c='black')
plt.xlabel(r'$\rm log_{10}(M_{\star}/M_{\odot})$')
plt.ylabel(r"$\tau_q$")
plt.legend()
plt.show()
#plt.title("bins of log_10 = 14.5 - 14.75 ")
#plt.savefig("plots_to_combine/all_bins_dex_split.pdf", format='pdf')

m_star = []
m_host = []

for i in range(7):

    p1 = np.where(data_t[i] == max(data_t[i]))[0]
    m_star.append(data_m[i][p1[0]])
    m_host.append(dex_range[i]+0.125)

from scipy.optimize import curve_fit

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B


A,B = curve_fit(f, m_host, np.log10(m_star))[0]
x = np.linspace(13,14.75)
fig = plt.figure(figsize=[7,5])
plt.scatter(m_host,np.log10(m_star), c="red", s=6)
plt.plot(x,A*x+B, c='black',label='{0:3f}*log_10(m_h)+{1:3f}'.format(A,B))
# del m_host[1]
# del m_star[1]
# A1,B1 = curve_fit(f,m_host, np.log10(m_star))[0]
# plt.plot(x,A1*x+B1, color = 'grey',label='{0:3f}*log_10(m_h)+{1:3f}, ignoring second entry'.format(A1,B1))
plt.ylabel(r'$\rm log_{10}(M_{\star}/M_{\odot})$')
plt.xlabel(r'$\rm log_{10}(M_{h}/M_{\odot})$')
plt.title(r'$ Location   \ of \  \ the \ peaks \ as \  a \ function \ of \ M_{h}$ ')
#plt.ylim(9.85,10.3)
#plt.legend()
# k =6
# p1 = np.where(data_t[k] == max(data_t[k]))[0]
# print(p1)
# print(np.log10(data_m[k][p1[0]]), m_host[k])
#plt.savefig("plots20.06.19/peak_as_a_f_m_h.pdf", format='pdf')
plt.show()
# print(p1)
# print(data_t[0][p1],data_m[0][p1])
#plt.show()
# print(host_gid[unique_hosts_index][h325])
# print(np.log10(M_h[unique_hosts_index][h325]))
# plt.hist(np.log10(M_h[unique_hosts_index]))
# plt.show()























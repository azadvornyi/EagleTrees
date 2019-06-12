import eagleSqlTools as sql
import numpy as np
from sys import argv
import matplotlib as mpl

mpl.use('TkAgg')
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
# ======================
# don't forget to switch to 3.3 rvir for the other simulations
# =====================
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

host_id = "21109760"

sim = "RefL0100N1504"
sat_id = 22

box_size = 100

host_gid, sat_num, M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504.txt", unpack=True,
                                                        usecols=(0, 1, 2, 3, 4, 5, 6))
# M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RecalL0100N0752.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))

n15 = np.where(M_h > 10 ** 15)
# print(len(M_h[n15]))
# print(max(M_h))
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
# print(m_s[keep_these])
a = set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
# small_hosts = np.where(M_h <= 1e14)
# small_hosts_set = set((small_hosts[0]))
# small_split = c.intersection(small_hosts_set)
# small_split = np.array(list(small_split))
small_split = c
#print(small_split)
#plt.rc('text', usetex=False)
with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f11:
    host_r_vir = pickle.load(f11)

with open("sim{0}/{1}/host_{1}".format(sim, host_id), 'rb') as f21:
    tree_host = pickle.load(f21)

with open("sim{0}/{1}/sat_{2}".format(sim, host_id, sat_id), 'rb') as f31:
    tree_sat = pickle.load(f31)

with open("sim{0}/{1}/sat_sub_{2}".format(sim, host_id, sat_id), 'rb') as f41:
    tree_sat_sub = pickle.load(f41)
print(len(tree_sat_sub['z']),len(tree_sat_sub['sub']))
for i in range(len(tree_sat_sub['sub'])):
    print(type(tree_sat_sub['sub'][i]))

sub_tree = tree_sat_sub['sub']
def times_Gyr(z):
    H0 = 67.77
    OmegaM = 0.307
    OmegaL = 0.693
    time_array = np.array([])
    for i in range(len(z)):
        t = (2 / (3 * H0 * np.sqrt(OmegaL))) * np.log(
            (np.sqrt(OmegaL * ((1 + z[i]) ** (-3))) + np.sqrt(OmegaL * ((1 + z[i]) ** (-3)) + OmegaM)) / np.sqrt(
                OmegaM))
        time_array = np.append(time_array, t * 1000)
    return time_array


def flag(scale_factor, box_, sat_pos, host_pos):
    coord = host_pos - sat_pos
    coord = coord * scale_factor
    half_box_ = box_ / 2
    if coord < -half_box_:
        coord = coord + box_
    elif coord > half_box_:
        coord = coord - box_
    return coord


# calculating distance between the host and the satellite
def moving_to_origin(host, sat, box, virial, sub):

    distances = np.array([])
    time_ = np.array([])
    perilist = np.array([])
    halfbox = box / 2
    if len(host) > len(sat):
        for i in reversed(range(len(sat))):
            for j in reversed(range(len(host))):

                if sat['z'][i] == host['z'][j]:
                    a = 1 / (1 + host['z'][j])
                    scaled_box = box * a
                    x = flag(a, scaled_box, sat['copx'][i], host['copx'][j])
                    y = flag(a, scaled_box, sat['copy'][i], host['copy'][j])
                    z = flag(a, scaled_box, sat['copz'][i], host['copz'][j])
                    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                    # virial = np.append(virial, r_v['z'][j])
                    distances = np.append(distances, dist)
                    time_ = np.append(time_, sat['z'][i])
    else:
        for i in reversed(range(len(host))):
            for j in reversed(range(len(sat))):

                if host['z'][i] == sat['z'][j]:
                    a = 1 / (1 + sat['z'][j])
                    scaled_halfbox = halfbox * a
                    scaled_box = box * a
                    x = sat['copx'][j] - host['copx'][i]
                    x = x * a
                    if x < -scaled_halfbox:
                        x = x + scaled_box
                    elif x > scaled_halfbox:
                        x = x - scaled_box

                    y = sat['copy'][j] - host['copy'][i]
                    y = y * a
                    if y < -scaled_halfbox:
                        y = y + scaled_box
                    elif y > scaled_halfbox:
                        y = y - scaled_box

                    z = sat['copz'][j] - host['copz'][i]
                    z = z * a
                    if z < -scaled_halfbox:
                        z = z + scaled_box
                    elif z > scaled_halfbox:
                        z = z - scaled_box

                    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)

                    # virial = np.append(virial, r_v['z'][i])
                    distances = np.append(distances, dist)
                    time_ = np.append(time_, host['z'][i])

    vir_time = virial['z'][::-1]
    vir_r = virial['r_vir'][::-1]  * 0.0015
    vir_r = vir_r*(1 + vir_time)
    sub = sub[::-1]
    sub1 = np.zeros(30)
    # host_mass = host_mass[::-1]
    # sat_mass = sat_mass[::-1]
    host_mass = np.zeros(30)
    sat_mass = np.zeros(30)
    distances = distances * (1 + time_)
    distances = distances
    # if sub is not None:

    if len(sat) > len(host):
        #print('top1')
        for i in range(len(vir_time)):
            #print('top2')
            for j in range(len(time_)):
                #print('top3')
                if vir_time[i] == time_[j]:
                    #print('top4')
                    if distances[j] < vir_r[i]:
                        #print('top5')


                        return host_mass[i],sat_mass[j], vir_time[i],sub[j],distances[j],distances,vir_r, vir_time,time_
                    # else:
                    #     #print('top6')
                    #     return -1, -1, -1, -1, -1, [-1], [-1], [-1], [-1]
    else:
        #print('bot1')
        for i in range(len(time_)):
            #print('bot2')
            for j in range(len(vir_time)):
                #print('bot3')
                if vir_time[j] == time_[i]:
                    #print('bot4')
                    if distances[i] < vir_r[j]:
                        #print('bot5')

                        return host_mass[i],sat_mass[i], vir_time[j],sub[i],distances[i],distances,vir_r,vir_time,time_

                    # else:
                    #     #print('bot6')
                    #     return -1, -1, -1, -1, -1, [-1], [-1], [-1], [-1]

    return -1, -1, -1, -1, -1, [-1], [-1], [-1], [-1]



# def passfunc(sub_T):
#     K =1
#     return K
# k = passfunc(sub_tree)


def a(z):
    return 1 / (1 + z)

print(type(np.array(None)))

sub_tree = np.array(sub_tree)
host_m, sat_m, zed,sub_n,d,rad,vir_rad, vir_time,time= moving_to_origin(tree_host, tree_sat, box_size, host_r_vir,
                                                                        sub_tree)
print(sub_n)
print(tree_sat_sub['sub'])
plt.plot(times_Gyr(vir_time),vir_rad,label = r'$\rm 1.5R_{vir}$')
plt.plot(times_Gyr(time),rad)
plt.scatter(times_Gyr([zed]),d, label = r'$ \rm  First \ timestep \ after \ crossing \ 1.5R_{vir}$')
plt.xlabel(r"$\rm t \ [Gyr]$")
plt.ylabel(r"$\rm R \ [cMpc]$")
plt.legend()
plt.savefig('first_time_count_after_crossing.pdf', format='pdf', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.show()
print(host_m, sat_m, zed,sub_n)
print(tree_host['mdm'],"mdm")
print(tree_sat['ms'],"ms")

print(tree_host['z'],'z_host')
print(tree_sat['z'],'z_sat')
#
sys.exit()


for i in small_split:

    with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
        host_r_vir = pickle.load(f1)

    with open("sim{0}/{1}/host_{1}".format(sim, int(host_gid[i])), 'rb') as f2:
        tree_host = pickle.load(f2)

    with open("sim{0}/{1}/sat_{2}".format(sim, int(host_gid[i]), int(sat_num[i])), 'rb') as f3:
        tree_sat = pickle.load(f3)
    with open("sim{0}/{1}/sat_sub_{2}".format(sim, int(host_id), int(sat_id)), 'rb') as f4:
        tree_sat_sub = pickle.load(f4)
    print(host_gid[i], sat_num[i])
    host_m, sat_m, zed, sub_n, d, rad,vir_rad, vir_time,time = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir, tree_sat_sub['sub'])
    print(sub_n)
    if sub_n != -1:
        f = open("the_ultimate_data_15_{0}".format(sim), "a")
        f.write(
            "{0} {1} {2} {3} {4} {5} {6} {7} {8} \n".format(int(host_gid[i]),int(sat_num[i]), tree_host['mdm'][0], max(tree_sat['ms']), t_i[i], t_q[i],
                                t_q[i]-t_i[i], sfr[i], sub_n))
        f.close()






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

host_id = "19054212"


sim = "RefL0100N1504"
sat_id = 7



box_size = 100


host_gid,sat_num,M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504_with_v.txt", unpack=True, usecols=(0,1,2, 3, 4, 5, 6))
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
# print(m_s[keep_these])
a = set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
small_hosts = np.where(M_h <= 1e14)
small_hosts_set = set((small_hosts[0]))
small_split = c.intersection(small_hosts_set)
small_split = np.array(list(small_split))
print(small_split)

with open("sim{0}_v/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as g1:
    host_r_vir = pickle.load(g1)

with open("sim{0}_v/{1}/host_{1}".format(sim, host_id), 'rb') as g2:
    tree_host = pickle.load(g2)

with open("sim{0}_v/{1}/sat_{2}".format(sim, host_id, sat_id), 'rb') as g3:
    tree_sat = pickle.load(g3)

with open("sim{0}_v/{1}/sat_vel_{2}".format(sim, host_id, sat_id), 'rb') as g4:
    tree_sat_v = pickle.load(g4)



#print(host_r_vir['v_z'])
#print(tree_sat['z'])
print(tree_host['mdm'])



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


def flag(scale_factor,box_,sat_pos,host_pos):
    coord = host_pos - sat_pos
    coord = coord * scale_factor
    half_box_ = box_/2
    if coord < -half_box_:
        coord = coord + box_
    elif coord > half_box_:
        coord = coord - box_
    return coord


# calculating distance between the host and the satellite
def moving_to_origin(host, sat, box, virial, sat_v):
    distances = np.array([])
    time_ = np.array([])
    perilist = np.array([])
    vir_holder = []
    halfbox = box / 2
    rel_v = []
    host_mass = []
    if len(host) > len(sat):
        for i in reversed(range(len(sat))):
            for j in reversed(range(len(host))):

                if sat['z'][i] == host['z'][j]:
                    a = 1 / (1 + host['z'][j])
                    scaled_box = box * a
                    x = flag(a,scaled_box,sat['copx'][i],host['copx'][j])
                    y = flag(a, scaled_box, sat['copy'][i], host['copy'][j])
                    z = flag(a, scaled_box, sat['copz'][i], host['copz'][j])
                    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                    v_s =   (virial['v_x'][j] - sat_v['v_x'][i]) ** 2 + (virial['v_y'][j] - sat_v['v_y'][i]) ** 2 + (virial['v_z'][j] - sat_v['v_z'][i]) ** 2
                    rel_v.append(v_s)
                    vir_holder.append(virial['r_vir'][j])
                    host_mass.append(host['mdm'][j])
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

                    v_s = (virial['v_x'][i] - sat_v['v_x'][j]) ** 2 + (virial['v_y'][i] - sat_v['v_y'][j]) ** 2 + (
                                virial['v_z'][i] - sat_v['v_z'][j]) ** 2
                    rel_v.append(v_s)
                    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                    vir_holder.append(virial['r_vir'][i])
                    host_mass.append(host['mdm'][i])
                    # virial = np.append(virial, r_v['z'][i])
                    distances = np.append(distances, dist)
                    time_ = np.append(time_, host['z'][i])

    # vir_time = virial['z'][::-1]
    # vir_r = virial['r_vir'][::-1]*0.0033
    print(host_mass)
    rel_v = np.array(rel_v)
    print(rel_v)
    vir_r = np.array(vir_holder)*0.0033
    host_mass = np.array(host_mass)
    # print(len(distances))
    # print(len(vir_r))
    distances = distances * (1 + time_)
    distances = distances
    time_ = time_
    time_ = times_Gyr(time_)
    peak, _ = sc.signal.find_peaks(-distances)
    for i in peak:
        if distances[i]<vir_r[i]:
            return distances, time_, i, time_[i], vir_r, rel_v, host_mass

    return distances, time_, -1, -1, -1, -1, -1



# radius, time_z, id, timeatid, virial_rad, v_rel, host_m = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir, tree_sat_v)
# #print(max(v_rel))
# peak, _ = sc.signal.find_peaks(-radius)
# print(peak)
# print(v_rel)
# print(type(v_rel))
# plt.plot(time_z,radius/max(radius))
# plt.plot(time_z, v_rel/max(v_rel))
# plt.plot(time_z, host_m/max(host_m))
# #plt.plot(time_z, virial_rad)
# plt.scatter(time_z[id], radius[id])
# plt.show()


for i in small_split:

    with open("sim{0}_v/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
        host_r_vir = pickle.load(f1)

    with open("sim{0}_v/{1}/host_{1}".format(sim, int(host_gid[i])), 'rb') as f2:
        tree_host = pickle.load(f2)

    with open("sim{0}_v/{1}/sat_{2}".format(sim, int(host_gid[i]), int(sat_num[i])), 'rb') as f3:
        tree_sat = pickle.load(f3)
    with open("sim{0}_v/{1}/sat_vel_{2}".format(sim, int(host_gid[i]), int(sat_num[i])), 'rb') as f4:
        tree_sat_v = pickle.load(f4)
    print(host_gid[i], sat_num[i])
    radius, time_z , pericenter_id, timeatid, virial_rad, v_rel, host_m = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir,tree_sat_v)
    # radius = np.array(radius)
    # for local minima
    #minima = argrelextrema(radius, np.less)
    #lessthan2mpc = np.where(radius<2)
    #firstset = set((minima[0]))
    #print(firstset)
    #secondset = set((lessthan2mpc[0]))
    # print(np.array(firstset.intersection(secondset)))
    #pericenter_id = (np.array(list(firstset.intersection(secondset))))
    #pericenter_id = minima[0]

    if (pericenter_id != -1)  and (radius[pericenter_id] != radius[-1]) and host_m[pericenter_id] != 0:


        rad = radius[pericenter_id]
        v_peri = v_rel[pericenter_id]
        mhost = host_m[pericenter_id]
        f = open("Rperi_1".format(sim), "a")
        f.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}\n".format(int(host_gid[i]), int(sat_num[i]),rad,M_h[i], m_s[i], t_i[i],
                                                                        t_q[i], sfr[i], np.max(tree_sat['ms']), v_peri, mhost))
        f.close()






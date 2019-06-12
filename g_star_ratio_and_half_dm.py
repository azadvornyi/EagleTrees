import numpy as np
from sys import argv
import matplotlib as mpl

mpl.use('TkAgg')
import matplotlib.pyplot as plt
import time
from scipy import interpolate
from sys import argv
import os
import pickle
import sys




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

box_size = 100
host_id = "20414028"


sim = "RefL0100N1504"
sat_id = 10

with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
    host_r_vir = pickle.load(f1)

with open("sim{0}/{1}/host_{1}".format(sim, host_id), 'rb') as f2:
    tree_host = pickle.load(f2)

with open("sim{0}/{1}/sat_{2}".format(sim, host_id, sat_id), 'rb') as f3:
    tree_sat = pickle.load(f3)


with open("sim{0}/{1}/sat_sub_{2}".format(sim, host_id, sat_id), 'rb') as f4:
    tree_sat_sub = pickle.load(f4)

print(tree_sat_sub['sub'])



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


def moving_to_origin(host, sat, box, r):
    distances = np.array([])
    time_ = np.array([])
    first_approach = np.array([])
    t_infall = np.array([])
    r_vir = np.array([])

    halfbox = box / 2
    if len(host) > len(sat):
        for i in reversed(range(len(sat))):
            for j in reversed(range(len(host))):

                if sat['z'][i] == host['z'][j]:
                    a = 1 / (1 + host['z'][j])
                    scaled_halfbox = halfbox * a
                    scaled_box = box * a
                    x = host['copx'][j] - sat['copx'][i]
                    x = x * a
                    if x < -scaled_halfbox:
                        x = x + scaled_box
                    elif x > scaled_halfbox:
                        x = x - scaled_box

                    y = host['copy'][j] - sat['copy'][i]
                    y = y * a
                    if y < -scaled_halfbox:
                        y = y + scaled_box
                    elif y > scaled_halfbox:
                        y = y - scaled_box

                    z = host['copz'][j] - sat['copz'][i]
                    z = z * a
                    if z < -scaled_halfbox:
                        z = z + scaled_box
                    elif z > scaled_halfbox:
                        z = z - scaled_box

                    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)

                    r_vir = np.append(r_vir, r[i])
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

                    r_vir = np.append(r_vir, r[i])
                    distances = np.append(distances, dist)
                    time_ = np.append(time_, host['z'][i])


    distances = list(reversed(distances))
    time_ = list(reversed(time_))
    time_ = times_Gyr(time_)

    return distances, time_


# print(tree_sat['ssfr'],'this is sat')



def t_infall_and_quench(r_vir, dist, time_dist, ssfr, gen_time, sat_time,mstar,mgas):

    """

    :param r_vir: 1
    :param dist: 2
    :param time_dist: 3
    :param ssfr: 5
    :param gen_time:
    :param sat_time:
    :return:
    """
    r_vir_func = interpolate.interp1d(gen_time, r_vir)
    dist_func = interpolate.interp1d(time_dist, dist)
    mgs_func = interpolate.interp1d(sat_time, mgas)
    mstar_func = interpolate.interp1d(sat_time,mstar)


    precise_time_d = np.linspace(min(time_dist), max(time_dist), 150)
    precise_time = np.linspace(min(gen_time), max(gen_time), 150)
    precise_time_sat = np.linspace(min(sat_time), max(sat_time), 150)

    interp_r = r_vir_func(precise_time)
    interp_d = dist_func(precise_time_d)
    interp_g = mgs_func(precise_time_sat)
    interp_m = mstar_func(precise_time_sat)

    first_approach = -1
    t_infall = -1
    t_quench = -1
    g_infall = -1
    m_infall = -1
    for counter in range(len(interp_r)):
        if interp_d[counter] < 0.0033 * interp_r[counter]:
            first_approach = interp_d[counter - 1]
            t_infall = precise_time_d[counter - 1]
            g_infall = interp_g[counter -1]
            m_infall = interp_m[counter -1]
            break

    ssfr_func = interpolate.interp1d(sat_time, ssfr)
    interp_ssfr = ssfr_func(precise_time_sat)
    for i in reversed(range(len(interp_ssfr))):
        if interp_ssfr[i] * 1e9 > 1e-2:
            t_quench = precise_time[i]
            break

    t_iq = t_quench - t_infall
    if first_approach == -1:
        return 0, 0, 0, 0

    return t_iq, t_infall, t_quench, first_approach, g_infall, m_infall, g_infall/m_infall


radius, time_z = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir['r_vir'])


tiq, t_infall, t_quench, first_approach_r, g_i, m_i, gm_i = t_infall_and_quench(host_r_vir['r_vir'], radius, time_z, tree_sat['ssfr'],
                                                                times_Gyr(tree_host['z']), times_Gyr((tree_sat['z'])),
                                                                tree_sat['ms'],tree_sat['m_g'])


print(times_Gyr(tree_sat['z']))
print(t_infall, g_i, m_i)
print(tree_sat['ms'])

sys.exit()


for i in small_split:

    with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
        host_r_vir = pickle.load(f1)

    with open("sim{0}/{1}/host_{1}".format(sim, int(host_gid[i])), 'rb') as f2:
        tree_host = pickle.load(f2)

    with open("sim{0}/{1}/sat_{2}".format(sim, int(host_gid[i]), int(sat_num[i])), 'rb') as f3:
        tree_sat = pickle.load(f3)
    #print(host_gid[i], sat_num[i])



    if (pericenter_id != -1)  and (radius[pericenter_id] != radius[-1]):


        rad = radius[pericenter_id]

        f = open("Rperi".format(sim), "a")
        f.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(int(host_gid[i]), int(sat_num[i]),rad,M_h[i], m_s[i], t_i[i], t_q[i], sfr[i], np.max(tree_sat['ms'])))
        f.close()
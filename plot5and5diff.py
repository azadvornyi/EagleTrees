
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
M_h, m_s, t_i, t_q, sfr = np.loadtxt("plot_info_for_short_and_long_sats_z2.txt", unpack=True, usecols=(2, 3, 4, 5, 6))



sat_flag = 1
with open("inspection/{0}/host_r_vir_{0}".format(sat_flag), 'rb') as g1:
    host_r_vir = pickle.load(g1)

with open("inspection/{0}/host_{0}".format(sat_flag), 'rb') as g2:
    tree_host = pickle.load(g2)

with open("inspection/{0}/sat_{0}".format(sat_flag), 'rb') as g3:
    tree_sat = pickle.load(g3)

with open("inspection/{0}/property1_sat_{0}".format(sat_flag), 'rb') as g4:
    p_sat_1 = pickle.load(g4)

with open("inspection/{0}/property2_sat_{0}".format(sat_flag), 'rb') as g5:
    p_sat_2 = pickle.load(g5)


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
def moving_to_origin(host, sat, box, virial):
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


                    dist = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                    vir_holder.append(virial['r_vir'][i])
                    host_mass.append(host['mdm'][i])
                    # virial = np.append(virial, r_v['z'][i])
                    distances = np.append(distances, dist)
                    time_ = np.append(time_, host['z'][i])

    # vir_time = virial['z'][::-1]

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
    return distances, time_, vir_r



box_size = 100

radius, time_z, virial_ = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir)
# #print(max(v_rel))
# peak, _ = sc.signal.find_peaks(-radius)
# print(peak)
# print(v_rel)
# print(type(v_rel))
# plt.plot(time_z,radius)
# plt.plot(time_z, virial_)
# plt.plot(time_z, v_rel/max(v_rel))
# plt.plot(time_z, host_m/max(host_m))
# #plt.plot(time_z, virial_rad)
# plt.scatter(time_z[id], radius[id])
#plt.show()

# sys.exit()
# for i in small_split:
#
#
#
#     if (pericenter_id != -1)  and (radius[pericenter_id] != radius[-1]) and host_m[pericenter_id] != 0:
#
#
#         rad = radius[pericenter_id]
#         v_peri = v_rel[pericenter_id]
#         mhost = host_m[pericenter_id]
#         f = open("Rperi_1".format(sim), "a")
#         f.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}\n".format(int(host_gid[i]), int(sat_num[i]),rad,M_h[i], m_s[i], t_i[i],
#                                                                         t_q[i], sfr[i], np.max(tree_sat['ms']), v_peri, mhost))
#         f.close()











# plot it
fig = plt.figure(figsize=(14, 12))
for i in np.arange(1,20):
    ins_1 = '1'
    ins_2 = '2'
    ins = ins_1
    sat_flag = i
    try:
        with open("inspection_z{1}/{0}/host_r_vir_{0}".format(sat_flag,ins), 'rb') as f1:
            host_r_vir = pickle.load(f1)

        with open("inspection_z{1}/{0}/host_{0}".format(sat_flag,ins), 'rb') as f2:
            tree_host = pickle.load(f2)

        with open("inspection_z{1}/{0}/sat_{0}".format(sat_flag,ins), 'rb') as f3:
            tree_sat = pickle.load(f3)

        with open("inspection_z{1}/{0}/property1_sat_{0}".format(sat_flag,ins), 'rb') as f4:
            p_sat_1 = pickle.load(f4)

        with open("inspection_z{1}/{0}/property2_sat_{0}".format(sat_flag,ins), 'rb') as f5:
            p_sat_2 = pickle.load(f5)
        alpha = 0.65
        if i <12:  #12 for z =2
            print(t_q[i-1]-t_i[i-1])
            color = 'red'
        else:
            color = 'blue'
            print(t_q[i-1] - t_i[i-1])
        radius, time_z, virial_ = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir)
        time_y = times_Gyr(tree_sat['z'])
        plt.subplot(431) # ssfr
        plt.semilogy(time_y-t_i[i-1],tree_sat['ssfr'], c= color,  alpha =alpha)
        plt.ylabel("ssfr")

        plt.subplot(432) # d_m
        plt.plot(time_y-t_i[i-1],tree_sat['mdm'], c= color, alpha =alpha)
        plt.ylabel("dark matter mass")
        plt.yscale('log')


        plt.subplot(433) # M_star
        plt.plot(time_y-t_i[i-1],tree_sat['ms'], c= color, alpha =alpha)
        plt.ylabel(r"$M_{star}$")
        plt.yscale('log')


        plt.subplot(434) # m_gas
        plt.plot(time_y-t_i[i-1],p_sat_1['m_g'], c= color, alpha =alpha)
        plt.ylabel(r"$M_{gas}$")
        plt.yscale('log')


        plt.subplot(435) # bhmar
        plt.plot(time_y-t_i[i-1],p_sat_1['bhmar'], c= color, alpha =alpha)
        plt.ylabel("Black hole mass accretion rate")
        plt.yscale('log')

        plt.subplot(436) # half_dm
        plt.plot(time_y-t_i[i-1],p_sat_1['half_dm'], c= color, alpha =alpha)
        plt.ylabel(r"$half mass R_{dm} $")


        plt.subplot(437) #half_stars
        plt.plot(time_y-t_i[i-1],p_sat_2['half_star'], c= color, alpha =alpha)
        plt.ylabel(r"$half mass R_{star}$")


        plt.subplot(438) #half_gas
        plt.plot(time_y-t_i[i-1],p_sat_2['half_gas'], c= color, alpha =alpha)
        plt.ylabel(r"$half mass R_{gas}$")
        plt.yscale('log')


        plt.subplot(439) #steellar velocity disp
        plt.plot(time_y-t_i[i-1],p_sat_2['svd'], c= color, alpha =alpha)
        plt.ylabel("stellar v dispersion")


        plt.subplot(4,3,12) #thermal energy
        plt.plot(time_z-t_i[i-1],radius, c= color, alpha =alpha)
        #plt.plot(time_z,virial_)
        plt.ylabel("distances ")
        plt.subplot(4,3,10)

        plt.text(0.5,0.8, r'$ z\approx 2, \  5.5<t_{i}/Gyr<6.5  $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
        plt.text(0.5, 0.3, r'$ \tau_{i} \approx 3.5 \ and  \ 8 \ Gyr  $',
                 rotation=0,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')
        plt.text(0.5, 0.5, r'$10^{9}<M_{sat}/M_{\odot}<10^{10} , \ \  2 \times 10^{13}< M_{host}/M_{\odot}<6 \times 10^{13} $',
                 rotation=0,
                 horizontalalignment='center',
                 verticalalignment='top',
                 multialignment='center')
        plt.tick_params(labelbottom=False, labelleft=False)
        plt.subplot(4,3,11) #thermal energy
        plt.plot(time_y-t_i[i-1],p_sat_1['ste'], c= color, alpha = alpha)
        plt.ylabel("specific thermal energy")
        plt.yscale('log')

        plt.tight_layout()
    except FileNotFoundError:
        pass

plt.savefig('plots20.06.19/10in1_z1.pdf', format='pdf', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.show()
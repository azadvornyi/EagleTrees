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
import pickle
import sys

# ======================
# don't forget to switch to 3.3 rvir for the other simulations
# =====================

host_id = "13995967"
sim = "RefL0100N1504"
sat_id = 2
box_size=100
with open("sim{0}_old/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
    host_r_vir = pickle.load(f1)

with open("sim{0}_old/{1}/host_{1}".format(sim, host_id), 'rb') as f2:
    tree_host = pickle.load(f2)

with open("sim{0}_old/{1}/sat_{2}".format(sim, host_id,sat_id), 'rb') as f3:
    tree_sat = pickle.load(f3)

def times_Gyr(z):
    H0 = 67.77
    OmegaM = 0.307
    OmegaL = 0.693
    time_array = np.array([])
    for i in range(len(z)):
        t = (2 / (3 * H0 * np.sqrt(OmegaL))) * np.log(
            (np.sqrt(OmegaL * ((1 + z[i]) ** (-3))) + np.sqrt(OmegaL * ((1 + z[i]) ** (-3)) + OmegaM)) / np.sqrt(
                OmegaM))
        time_array = np.append(time_array,t*1000)
    return time_array


# calculating distance between the host and the satellite
def moving_to_origin(host, sat, box, r):
    distances = np.array([])
    time_ = np.array([])
    first_approach = np.array([])
    t_infall = np.array([])
    r_vir = np.array([])

    halfbox = box/2
    if len(host) > len(sat):
        for i in reversed(range(len(sat))):
            for j in reversed(range(len(host))):

                if sat['z'][i] == host['z'][j]:
                    a = 1/(1+host['z'][j])
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

                    dist = np.sqrt(x**2 + y**2 + z**2)

                    r_vir = np.append(r_vir, r[i])
                    distances = np.append(distances,dist)
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


    for counter in range(len(distances)):
        if distances[counter] < 0.0025 * r_vir[counter]:
            first_approach = np.append(first_approach,distances[counter])
            t_infall = np.append(t_infall, time_[counter])
    distances = list(reversed(distances))
    time_ = list(reversed(time_))
    time_ = times_Gyr(time_)

    return distances, time_


#print(tree_sat['ssfr'],'this is sat')

radius, time_z = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir['r_vir'])

#print(time.time() - t, 'it took this many seconds')


# calculating t_quench - t_infall
def t_infall_and_quench(r_vir, dist, time_dist, ssfr, gen_time, sat_time):
    r_vir_func = interpolate.interp1d(gen_time, r_vir )
    dist_func = interpolate.interp1d(time_dist, dist )

    precise_time_d = np.linspace(min(time_dist),max(time_dist), 150)
    precise_time = np.linspace(min(gen_time), max(gen_time), 150)
    precise_time_sat = np.linspace(min(sat_time), max(sat_time), 150)

    interp_r = r_vir_func(precise_time)
    interp_d = dist_func(precise_time_d)

    first_approach = -1
    t_infall = -1
    t_quench = -1
    for counter in range(len(interp_r)):
        if interp_d[counter] < 0.0025 * interp_r[counter]:
            first_approach = interp_d[counter-1]
            t_infall = precise_time_d[counter-1]
            break


    ssfr_func = interpolate.interp1d(sat_time,ssfr)
    interp_ssfr = ssfr_func(precise_time_sat)
    for i in reversed(range(len(interp_ssfr))):
        if interp_ssfr[i]*1e9>1e-2:
            t_quench = precise_time[i]
            break
        else:
            pass

    t_iq = t_quench - t_infall
    if first_approach == -1:
        return 0, 0, 0, 0
    else:
        return t_iq, t_infall, t_quench, first_approach


tiq, t_infall, t_quench, first_approach_r = t_infall_and_quench(host_r_vir['r_vir'], radius,  time_z, tree_sat['ssfr'],
                                                                times_Gyr(tree_host['z']),times_Gyr((tree_sat['z'])))



# plt.savefig('dist.pdf', format ='pdf')
# plt.clf()
# ========right side ssfr plots
fig = plt.figure(figsize=[12,9])
plt.subplot(321)
plt.plot(time_z, radius, c = 'blue')
plt.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0025,c = 'green')
plt.hlines(first_approach_r,0,max(time_z),colors='orange',linestyles='--')
plt.vlines(t_infall,-0.1, max(radius)+0.1, colors='orange',linestyles='--')
plt.tick_params(direction='in', length=4, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
plt.autoscale(enable=False, axis='y', tight=None)
plt.twiny()
#reversedz =np.fliplr(host_r_vir['z'])
reversed_arr = host_r_vir['z'][::-1]
print(reversed_arr)
#print(host_r_vir['z'])
plt.plot(host_r_vir['z'],host_r_vir['r_vir'][::-1]*0.0025)
#plt.grid()
plt.subplot(323)
plt.plot(time_z, radius, c = 'blue', label="radial distance of a satellite")
plt.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0025,c = 'green', label = "3.3R_{vir} of a host [cMpc]")
plt.ylabel('radius, [cMpc]')
plt.legend()
plt.hlines(first_approach_r,0,max(time_z),colors='orange',linestyles='--')
plt.vlines(t_infall,-0.1, max(radius)+0.1, colors='orange',linestyles='--')
plt.tick_params(direction='in', length=4, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
#plt.grid()

plt.subplot(325)
plt.plot(time_z, radius, c = 'blue')
plt.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0025,c = 'green')

plt.hlines(first_approach_r,0,max(time_z),colors='orange',linestyles='--')
plt.vlines(t_infall,-0.1, max(radius) + 0.1, colors='orange',linestyles='--')
plt.tick_params(direction='in', length=4, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
#plt.grid()
plt.xlabel('time since BigBang [Gyr]')
plt.subplots_adjust(hspace = 0)
# plt.ylabel("t_q [Gyr]")
# plt.xlabel("t_i [Gyr]")
#
# plt.subplot(323)
# plt.hist(t_i[c],bins=28)
# plt.xlabel("t_i [Gyr]")
# plt.ylabel("N")
# #================
# plt.subplot(325)
# plt.scatter(m_s[c], t_q[c]-t_i[c],s= 2, c= "blue")
# plt.xscale('log')
# plt.ylabel("t_q - t_i [Gyr]")
# plt.xlabel("M_sat [M_sun]")
# # =================left side distance plots
# plt.subplot(322)
# plt.hist(np.log(m_s[c]), bins=28)
# plt.xscale('log')
# plt.xlabel("M_sat [M_sun]")
# plt.ylabel("N")
# #==============
# plt.subplot(324)
# plt.scatter(t_i[c], sfr[c]*1e9, s= 2, c= "blue")
# plt.ylabel("ssfr(z=0)  Gyr^-1 ")
# plt.xlabel("t_i [Gyr]")
# # plt.yscale('log')
# # plt.ylim(1e-9,1)
# plt.tight_layout()
#
# plt.subplot(326)
# plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=1,
#            label='binned statistic of data')
# plt.ylabel("Quenching timescale [Gyr]")
# #plt.plot((binnumber - 0.5) * bin_width, x_pdf, 'g.', alpha=0.5)
# plt.xlabel("M_sat [M_sun]")

# Do the plot code
plt.savefig('ssfr_dist.pdf', format='pdf', dpi=1200, pad_inches = 0.1, bbox_inches = 'tight')
plt.show()
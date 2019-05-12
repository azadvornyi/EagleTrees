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
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.cosmology import z_at_value
# ======================
# don't forget to switch to 3.3 rvir for the other simulations
# =====================

host_id = "13995967"
host_id2 = "13995967"
host_id3 = "13995967"
sim = "RefL0100N1504"
sat_id = 2
sat_id2 = 2
sat_id3 = 2
box_size=100
with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
    host_r_vir = pickle.load(f1)

with open("sim{0}/{1}/host_{1}".format(sim, host_id), 'rb') as f2:
    tree_host = pickle.load(f2)

with open("sim{0}/{1}/sat_{2}".format(sim, host_id,sat_id), 'rb') as f3:
    tree_sat = pickle.load(f3)


#
# with open("sim{0}_old/{1}/host_r_vir_{1}".format(sim, host_id2), 'rb') as f4:
#     host_r_vir2 = pickle.load(f4)
#
# with open("sim{0}_old/{1}/host_{1}".format(sim, host_id2), 'rb') as f5:
#     tree_host2 = pickle.load(f5)
#
# with open("sim{0}_old/{1}/sat_{2}".format(sim, host_id2,sat_id), 'rb') as f6:
#     tree_sat2 = pickle.load(f6)
#
#
#
# with open("sim{0}_old/{1}/host_r_vir_{1}".format(sim, host_id3), 'rb') as f7:
#     host_r_vir3 = pickle.load(f7)
#
# with open("sim{0}_old/{1}/host_{1}".format(sim, host_id3), 'rb') as f8:
#     tree_host3 = pickle.load(f8)
#
# with open("sim{0}_old/{1}/sat_{2}".format(sim, host_id,sat_id3), 'rb') as f9:
#     tree_sat3 = pickle.load(f9)

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
        if interp_d[counter] < 0.0033 * interp_r[counter]:
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



cosmo = FlatLambdaCDM(H0=67.77, Om0=0.307)


redshifts_to_display = np.array([0,0.5,1,2,4,100])
#ageticks = [z_at_value(cosmo.age, age) for age in ages]
ages = times_Gyr(redshifts_to_display)#
#first
font = {'family' : 'monospace',
        'size'   : '12'}
header = {'family' : 'monospace',
        'size'   : '15'}
plt.rc('font', **font)
tmax, tmin = 14.5, -0.5
fig = plt.figure(figsize=(14,9))
ax = fig.add_subplot(321)
#ax.plot(time_z, radius * (1+tree_sat['z'])) #comoving
#ax.set_xticks(ages.value)
ax2 = ax.twiny()

ax2.set_xticks(ages)

ax2.set_xticklabels(('0','0.5','1','2','4',r'$\infty$'))


ax2.set_xlim(tmin, tmax)
ax2.set_xlabel('Redshift z')
ax2.tick_params(direction='in', length=4, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)

ax.plot(time_z, radius, c = 'blue')

ax.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0033,c = 'green')
ax.hlines(first_approach_r,0,max(time_z),colors='black',linestyles='--')
ax.vlines(t_infall,-0.1, max(radius)+0.1, colors='black',linestyles='--')
ax.text(16.6, first_approach_r+0.055, r'$  \longleftarrow \mathtt{3.3R_{200crit}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax.text(t_infall+1.5, first_approach_r+0.4, r'$  \longleftarrow \mathtt{t_{infall}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
#ax.set_title('Orbital history')
ax.text(7, 2.5, r'$\mathtt{Orbital \ \ \ history}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center',
        fontdict=header)
ax.tick_params(direction='in', length=0, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
ax.set_xlim(tmin, tmax)

ax3 = fig.add_subplot(323)
ax3.plot(time_z, radius, c = 'blue', label=r'$\mathtt{R_{sep,sat}}$')
ax3.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0033,c = 'green',label=r'$\mathtt{R_{vir,h}}$')
ax3.hlines(first_approach_r,0,max(time_z),colors='black',linestyles='--')
ax3.vlines(t_infall,-0.1, max(radius)+0.1, colors='black',linestyles='--')
ax3.tick_params(direction='in', length=0, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
ax3.set_xlim(tmin, tmax)

ax3.text(16.6, first_approach_r+0.055, r'$  \longleftarrow \mathtt{3.3R_{200crit}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax3.text(t_infall+1.5, first_approach_r+0.4, r'$  \longleftarrow \mathtt{t_{infall}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')


plt.legend()
ax4 = fig.add_subplot(325)
ax4.plot(time_z, radius, c = 'blue')
ax4.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0033,c = 'green')
ax4.hlines(first_approach_r,0,max(time_z),colors='black',linestyles='--')
ax4.vlines(t_infall,-0.1, max(radius)+0.1, colors='black',linestyles='--')
ax4.tick_params(direction='in', length=4, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
ax4.set_xlim(tmin, tmax)

ax4.text(16.6, first_approach_r+0.055, r'$  \longleftarrow \mathtt{3.3R_{200crit}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax4.text(t_infall+1.5, first_approach_r+0.4, r'$  \longleftarrow \mathtt{t_{infall}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')




ax4.set_xlabel('time since Big Bang [Gyr]')
ax3.set_ylabel(r'$\mathtt{R_{host} \ \ [pMpc]}$')



ax = fig.add_subplot(322)
#ax.plot(time_z, radius * (1+tree_sat['z'])) #comoving
#ax.set_xticks(ages.value)
ax2 = ax.twiny()

ax2.set_xticks(ages)

ax2.set_xticklabels(('0','0.5','1','2','4',r'$\infty$'))


ax2.set_xlim(tmin, tmax)
ax2.set_xlabel('Redshift z')
ax2.tick_params(direction='in', length=4, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)

ax.semilogy(times_Gyr(tree_sat['z']), tree_sat['ssfr'], '-r')
ax.vlines(t_quench,min(tree_sat['ssfr']), max(tree_sat['ssfr']) ,colors='black',linestyles='--')
ax.hlines(tree_sat['ssfr'][0],min(time_z),max(time_z),colors='black',linestyles='--')
#ax.vlines(t_infall,-0.1, max(radius)+0.1, colors='orange',linestyles='--')
ax.text(12.5, 10**(-9.5), r'$  \mathtt{t_{quench}} \longrightarrow $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax.text(7, tree_sat['ssfr'][0]+2*10**(-11), r'$ \mathtt{SSFRcut}  $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax.text(7, 10**(-6.85), r'$\mathtt{SSFR \ \ \ history}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center',
        fontdict=header)
ax.tick_params(direction='in', length=0, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
ax.set_xlim(tmin, tmax)

ax3 = fig.add_subplot(324)
ax3.plot(time_z, radius, c = 'blue')
ax3.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0025,c = 'green')
ax3.hlines(first_approach_r,0,max(time_z),colors='black',linestyles='--')
ax3.vlines(t_infall,-0.1, max(radius)+0.1, colors='black',linestyles='--')
ax3.tick_params(direction='in', length=0, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
ax3.set_xlim(tmin, tmax)

ax3.text(16.6, first_approach_r+0.055, r'$  \longleftarrow \mathtt{3.3R_{200crit}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax3.text(t_infall+1.5, first_approach_r+0.4, r'$  \longleftarrow \mathtt{t_{infall}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')



ax4 = fig.add_subplot(326)
ax4.plot(time_z, radius, c = 'blue')
ax4.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0025,c = 'green')
ax4.hlines(first_approach_r,0,max(time_z),colors='black',linestyles='--')
ax4.vlines(t_infall,-0.1, max(radius)+0.1, colors='black',linestyles='--')
ax4.tick_params(direction='in', length=4, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
ax4.set_xlim(tmin, tmax)

ax4.text(16.6, first_approach_r+0.055, r'$  \longleftarrow \mathtt{3.3R_{200crit}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax4.text(t_infall+1.5, first_approach_r+0.4, r'$  \longleftarrow \mathtt{t_{infall}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')




ax4.set_xlabel('time since BigBang [Gyr]')
ax3.set_ylabel(r'$\mathtt{log(SSFR) [yr^{-1}]}$')
plt.subplots_adjust(hspace = 0, wspace = 0.48)




# Do the plot code
plt.savefig('ssfr_dist.pdf', format='pdf', dpi=1200, pad_inches = 0.1, bbox_inches = 'tight')
plt.show()
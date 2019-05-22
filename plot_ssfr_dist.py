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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.lines as mlines
import astropy.units as u
from astropy.cosmology import z_at_value
# ======================
# don't forget to switch to 3.3 rvir for the other simulations
# =====================
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

host_id = "13995967"
host_id2 = "18737995" # host_id2 = "19634929"
host_id3 = "14071649"

sim = "RefL0100N1504"
sat_id = 2
sat_id2 = 35 #sat_id2 = 41
sat_id3 = 2

box_size = 100

with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
    host_r_vir = pickle.load(f1)

with open("sim{0}/{1}/host_{1}".format(sim, host_id), 'rb') as f2:
    tree_host = pickle.load(f2)

with open("sim{0}/{1}/sat_{2}".format(sim, host_id, sat_id), 'rb') as f3:
    tree_sat = pickle.load(f3)

with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id2), 'rb') as f4:
    host_r_vir2 = pickle.load(f4)

with open("sim{0}/{1}/host_{1}".format(sim, host_id2), 'rb') as f5:
    tree_host2 = pickle.load(f5)

with open("sim{0}/{1}/sat_{2}".format(sim, host_id2, sat_id2), 'rb') as f6:
    tree_sat2 = pickle.load(f6)

with open("sim{0}/{1}/host_r_vir_{1}".format(sim, host_id3), 'rb') as f7:
    host_r_vir3 = pickle.load(f7)

with open("sim{0}/{1}/host_{1}".format(sim, host_id3), 'rb') as f8:
    tree_host3 = pickle.load(f8)

with open("sim{0}/{1}/sat_{2}".format(sim, host_id3, sat_id3), 'rb') as f9:
    tree_sat3 = pickle.load(f9)

print(tree_sat2['ssfr'])
print(tree_sat2['z'])


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

print(times_Gyr((tree_sat2['z'])))
# calculating distance between the host and the satellite
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


    distances = distances[::-1]
    time_ = time_[::-1]
    time_ = times_Gyr(time_)

    return distances, time_


# print(tree_sat['ssfr'],'this is sat')

radius, time_z = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir['r_vir'])
radius2, time_z2 = moving_to_origin(tree_host2, tree_sat2, box_size, host_r_vir2['r_vir'])
radius3, time_z3 = moving_to_origin(tree_host3, tree_sat3, box_size, host_r_vir3['r_vir'])


# print(time.time() - t, 'it took this many seconds')


# calculating t_quench - t_infall
def t_infall_and_quench(r_vir, dist, time_dist, ssfr, gen_time, sat_time):
    r_vir_func = interpolate.interp1d(gen_time, r_vir)
    dist_func = interpolate.interp1d(time_dist, dist)

    precise_time_d = np.linspace(min(time_dist), max(time_dist), 150)
    precise_time = np.linspace(min(gen_time), max(gen_time), 150)
    precise_time_sat = np.linspace(min(sat_time), max(sat_time), 150)

    interp_r = r_vir_func(precise_time)
    interp_d = dist_func(precise_time_d)

    first_approach = -1
    t_infall = -1
    t_quench = -1
    for counter in range(len(interp_r)):
        if interp_d[counter] < 0.0033 * interp_r[counter]:
            first_approach = interp_d[counter - 1]
            t_infall = precise_time_d[counter - 1]
            break

    ssfr_func = interpolate.interp1d(sat_time, ssfr)
    interp_ssfr = ssfr_func(precise_time_sat)
    for i in reversed(range(len(interp_ssfr))):
        if interp_ssfr[i] * 1e9 > 1e-2:
            t_quench = precise_time[i]
            break
        else:
            pass

    t_iq = t_quench - t_infall
    if first_approach == -1:
        return 0, 0, 0, 0
    else:
        return t_iq, t_infall, t_quench, first_approach


tiq, t_infall, t_quench, first_approach_r = t_infall_and_quench(host_r_vir['r_vir'], radius, time_z, tree_sat['ssfr'],
                                                                times_Gyr(tree_host['z']), times_Gyr((tree_sat['z'])))
tiq2, t_infall2, t_quench2, first_approach_r2 = t_infall_and_quench(host_r_vir2['r_vir'], radius2, time_z2,
                                                                    tree_sat2['ssfr'],
                                                                    times_Gyr(tree_host2['z']),
                                                                    times_Gyr((tree_sat2['z'])))
tiq3, t_infall3, t_quench3, first_approach_r3 = t_infall_and_quench(host_r_vir3['r_vir'], radius3, time_z3,
                                                                    tree_sat3['ssfr'],
                                                                    times_Gyr(tree_host3['z']),
                                                                    times_Gyr((tree_sat3['z'])))


def a(z):
    return 1 / (1 + z)


cosmo = FlatLambdaCDM(H0=67.77, Om0=0.307)

redshifts_to_display = np.array([0, 0.5, 1, 2, 4, 10000])
# ageticks = [z_at_value(cosmo.age, age) for age in ages]
ages = times_Gyr(redshifts_to_display)  #
# first
font = {'family': 'monospace',
        'size': '12'}
header = {'family': 'monospace',
          'size': '15'}
small_font = {'family': 'monospace',
              'size': '10'}
plt.rc('font', **font)

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(5)

tmax, tmin = 14.3, -0.05
fig = plt.figure(figsize=(16, 10))
# +++++++++=========-==-=-=-=--=--==-=-=-=-=-=--==-=-=-=-=-=-=-=-=-=-=-=--=-==-=--=-=-==-=-=-=-
ax = fig.add_subplot(321)
# ax.plot(time_z, radius * (1+tree_sat['z'])) #comoving
# ax.set_xticks(ages.value)
ax2 = ax.twiny()

ax2.set_xticks(ages)

ax2.set_xticklabels(('0', '0.5', '1', '2', '4', r'$\infty$'))

ax2.set_xlim(tmin, tmax)
ax2.set_xlabel('Redshift z')
ax2.tick_params(direction='in', length=6, width=1.1, colors='black',
                grid_color='black', grid_alpha=0.4)

ax.plot(time_z, radius / a(tree_sat['z']), c='blue', zorder=2)
ax.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir'] * 0.0033 / a(host_r_vir['z']), c='green', zorder=1)
ax.hlines(1.86, -1, 15, colors='black', linestyles='--', lw=1, zorder=10)
ax.vlines(t_infall, -0.1, 3.5, colors='black', linestyles='--', lw=1, zorder=10)
ax.text(15.3, 1.955, r'$  \longleftarrow \mathtt{R_{c}}$',
        rotation=0,
        horizontalalignment='center',
        verticalalignment='top',
        multialignment='center')
ax.text(t_infall + 1.2, 2.8, r'$  \longleftarrow \mathtt{t_{infall}}$',
        rotation=0,
        horizontalalignment='center',
        verticalalignment='top',
        multialignment='center')
# ax.set_title('Orbital history')
ax.text(7, 4.9, r'$\mathtt{Orbital \ \ \ history}$',
        rotation=0,
        horizontalalignment='center',
        verticalalignment='top',
        multialignment='center',
        fontdict=header)
ax.tick_params(direction='in', length=0, width=1, colors='black',
               grid_color='black', grid_alpha=0.4)
ax.set_xlim(tmin, tmax)
ax.set_ylim(0, 3.5)
minorLocatory = AutoMinorLocator()
minorLocatorx = AutoMinorLocator()
ax.xaxis.set_minor_locator(minorLocatorx)
ax.yaxis.set_minor_locator(minorLocatory)
ax.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
               grid_color='black', grid_alpha=0.4, labelbottom=False)
ax.tick_params(axis='x', which='minor', length=4, color='black', direction='in')
ax.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax.tick_params(axis='y', direction='in', length=6, width=1.1, which='major', colors='black',
               grid_color='black', grid_alpha=0.4)

ax3 = fig.add_subplot(323)
ax2 = ax3.twiny()

ax2.set_xticks(ages)
ax2.set_xticklabels((' ', ' ', ' ', ' ', ' ', ' '))

ax2.set_xlim(tmin, tmax)
ax2.tick_params(direction='in', length=6, width=1.1, colors='black',
                grid_color='black', grid_alpha=0.4)

ax3.plot(time_z2, radius2 / a(host_r_vir2['z']), c='blue')
ax3.plot(times_Gyr(host_r_vir2['z']), host_r_vir2['r_vir'] * 0.0033 / a(host_r_vir2['z']), c='green')
ax3.hlines(3.15, -1, 15, colors='black', linestyles='--', lw=1, zorder=4)
ax3.vlines(t_infall2, -1, 9, colors='black', linestyles='--', lw=1, zorder=4)
ax3.tick_params(direction='in', length=0, width=1, colors='black',
                grid_color='black', grid_alpha=0.4)
ax3.set_xlim(tmin, tmax)
ax3.set_ylim(0, 7.5)
ax3.text(15.3, 3.36, r'$  \longleftarrow \mathtt{R_{c}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax3.text(t_infall2 + 1.2, 5, r'$    \longleftarrow \mathtt{t_{infall}}  $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
minorLocatory = AutoMinorLocator()
minorLocatorx = AutoMinorLocator()
ax3.xaxis.set_minor_locator(minorLocatorx)
ax3.yaxis.set_minor_locator(minorLocatory)
ax3.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4, labelbottom=False)
ax3.tick_params(axis='x', which='minor', length=4, color='black', direction='in')
ax3.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax3.tick_params(axis='y', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4)
blue_line = mlines.Line2D([], [], color='green', label=r'$\mathtt{3.3R_{200c,host}}$')
green_line = mlines.Line2D([], [], color='blue', label=r'$\mathtt{R_{sep,sat}}$')

plt.legend(handles=[blue_line, green_line], frameon=False)
ax4 = fig.add_subplot(325)

ax2 = ax4.twiny()

ax2.set_xticks(ages)
ax2.set_xticklabels((' ', ' ', ' ', ' ', ' ', ' '))

ax2.set_xlim(tmin, tmax)
ax2.tick_params(direction='in', length=6, width=1.1, colors='black',
                grid_color='black', grid_alpha=0.4)
ax4.plot(time_z3, radius3 / a(tree_sat3['z']), c='blue', label=r'$\mathtt{R_{sep,sat}}$')
ax4.plot(times_Gyr(host_r_vir3['z']), host_r_vir3['r_vir'] * 0.0033 / a(host_r_vir3['z']), c='green',
         label=r'$\mathtt{3.3R_{200c,h}}$')
ax4.hlines(1.605, -1, 15, colors='black', linestyles='--', lw=1, zorder=5)
ax4.vlines(t_infall3, -0.1, 4, colors='black', linestyles='--', lw=1, zorder=5)
ax4.tick_params(direction='in', length=0, width=1, colors='black',
                grid_color='black', grid_alpha=0.4)
ax4.set_xlim(tmin, tmax)

ax4.text(15.3, 1.68, r'$  \longleftarrow \mathtt{R_{c}}$',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax4.text(t_infall3 + 1.2, 2.3, r'$  \longleftarrow \mathtt{t_{infall}} $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
minorLocatory = AutoMinorLocator()
minorLocatorx = AutoMinorLocator()
ax4.xaxis.set_minor_locator(minorLocatorx)
ax4.yaxis.set_minor_locator(minorLocatory)
ax4.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4)
ax4.tick_params(axis='x', which='minor', length=4, color='black', direction='in')
ax4.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax4.tick_params(axis='y', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4)

ax4.set_ylim(0, 2.9)
ax4.set_xlabel('Time since Big Bang [Gyr]')
ax3.set_ylabel(r'$\mathtt{R_{host} \ \ [cMpc]}$')

# +++++++++=========-==-=-=-=--=--==-=-=-=-=-=--==-=-=-=-=-=-=-=-=-=-=-=--=-==-=--=-=-==-=-=-=-
minorLocatory = AutoMinorLocator()
minorLocatorx = AutoMinorLocator()
ax = fig.add_subplot(322)
# ax.plot(time_z, radius * (1+tree_sat['z'])) #comoving
# ax.set_xticks(ages.value)
ax2 = ax.twiny()

ax2.set_xticks(ages)

ax2.set_xticklabels(('0', '0.5', '1', '2', '4', r'$\infty$'))

ax2.set_xlim(tmin, tmax)
ax2.set_xlabel('Redshift z')
ax2.tick_params(direction='in', length=4, width=1, colors='black',
                grid_color='black', grid_alpha=0.4)

ax.plot(times_Gyr(tree_sat['z']), (tree_sat['ssfr']), '-r')
ax.set_yscale('log')
ax.vlines(t_quench, min(tree_sat['ssfr']), max(tree_sat['ssfr'] + 1), colors='black', linestyles='--', lw=1, zorder=10)
ax.hlines(tree_sat['ssfr'][0], tmin, tmax, colors='black', linestyles='--', lw=1, zorder=10)
# ax.vlines(t_infall,-0.1, max(radius)+0.1, colors='orange',linestyles='--')
ax.text(t_quench - 1.25, 10 ** (-9.5), r'$  \mathtt{t_{quench}} \longrightarrow $',
        rotation=0,
        horizontalalignment='center',
        verticalalignment='top',
        multialignment='center')
ax.text(7, tree_sat['ssfr'][0] + 2 * 10 ** (-11), r'$ \mathtt{SSFRcut}  $',
        rotation=0,
        horizontalalignment='center',
        verticalalignment='top',
        multialignment='center')
ax.text(7, 10 ** (-5.7), r'$\mathtt{SSFR \ \ \ history}$',
        rotation=0,
        horizontalalignment='center',
        verticalalignment='top',
        multialignment='center',
        fontdict=header)
ax.xaxis.set_minor_locator(minorLocatorx)

ax.xaxis.set_minor_locator(minorLocatorx)

ax.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
               grid_color='black', grid_alpha=0.4, labelbottom=False)
ax.tick_params(which='minor', length=4, color='black', direction='in')
ax.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
               grid_color='black', grid_alpha=0.4)

ax.set_xlim(tmin, tmax)
ax2.set_ylim(5 * 10 ** (-12), 5 * 10 ** (-8))

ax3 = fig.add_subplot(324)
ax2 = ax3.twiny()

ax2.set_xticks(ages)
ax2.set_xticklabels((' ', ' ', ' ', ' ', ' ', ' '))

ax2.set_xlim(tmin, tmax)
ax2.tick_params(direction='in', length=6, width=1.1, colors='black',
                grid_color='black', grid_alpha=0.4)
ax3.set_ylim(1.5 * 10 ** (-12), 5 * 10 ** (-8))
ax3.plot(times_Gyr(tree_sat2['z']), (tree_sat2['ssfr']), '-r')
ax3.set_yscale('log')
ax3.vlines(t_quench2, min(tree_sat2['ssfr']), max(tree_sat2['ssfr'] + 0.5), colors='black', linestyles='--', lw=1,
           zorder=10)
ax3.hlines(tree_sat2['ssfr'][1] + 3 * 10 ** (-12), tmin, tmax, colors='black', linestyles='--', lw=1, zorder=10)
# ax.vlines(t_infall,-0.1, max(radius)+0.1, colors='orange',linestyles='--')
ax3.text(t_quench2 - 1.25, 10 ** (-9.5), r'$  \mathtt{t_{quench}} \longrightarrow $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax3.text(7, tree_sat2['ssfr'][0] + 1.77 * 10 ** (-11), r'$ \mathtt{SSFRcut}  $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')

red_line = mlines.Line2D([], [], color='red', label=r"$\mathtt{SSFR_{sat}} $")

plt.legend(handles=[red_line], frameon=False, loc=9)
# plt.legend(loc = 9,fancybox = True, frameon=False)
ax3.set_xlim(tmin, tmax)
ax3.xaxis.set_minor_locator(minorLocatorx)

ax3.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4, labelbottom=False)
ax3.tick_params(which='minor', length=4, color='black', direction='in')
ax3.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax3.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
                grid_color='black', grid_alpha=0.4)

minorLocatory = AutoMinorLocator()
minorLocatorx = AutoMinorLocator()
ax4 = fig.add_subplot(326)
ax2 = ax4.twiny()

ax2.set_xticks(ages)
ax2.set_xticklabels((' ', ' ', ' ', ' ', ' ', ' '))

ax2.set_xlim(tmin, tmax)
ax2.tick_params(direction='in', length=6, width=1.1, colors='black',
                grid_color='black', grid_alpha=0.4)
ax4.set_ylim(2 * 10 ** (-12), 5 * 10 ** (-8))
ax4.plot(times_Gyr(tree_sat3['z']), (tree_sat3['ssfr']), '-r')
ax4.set_yscale('log')
ax4.vlines(t_quench3, min(tree_sat3['ssfr']), max(tree_sat3['ssfr'] + 0.5), colors='black', linestyles='--', lw=1,
           zorder=10)
ax4.hlines(10 ** (-11), tmin, tmax, colors='black', linestyles='--', lw=1, zorder=10)
# ax.vlines(t_infall,-0.1, max(radius)+0.1, colors='orange',linestyles='--')
ax4.text(t_quench3 + 1.25, 2.8 * 10 ** (-11), r'$ \longleftarrow \mathtt{t_{quench}}  $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
ax4.text(5, 2.3 * 10 ** (-11), r'$ \mathtt{SSFRcut}  $',
         rotation=0,
         horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')

ax4.xaxis.set_minor_locator(minorLocatorx)

ax4.tick_params(axis='x', direction='in', length=6, width=1.1, which='major', colors='black',
                grid_color='black', grid_alpha=0.4)
ax4.tick_params(which='minor', length=4, color='black', direction='in')
ax4.tick_params(axis='y', which='minor', length=4, color='black', direction='in')
ax4.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
                grid_color='black', grid_alpha=0.4)

ax4.set_xlim(tmin, tmax)
# ax4.vlines(12,0,1)
# ax4.vlines(10,0,1)
# ax4.vlines(8,0,1)

# axins = inset_axes(ax4, width=1.855, height=1.1)
axins = inset_axes(ax4, width="31.42%", height="50%",
                   bbox_to_anchor=(.5456, .46, .885, .9),
                   bbox_transform=ax4.transAxes, loc=3)
axins.tick_params(direction='in', labelbottom=False)  #
axins.plot(times_Gyr(tree_sat3['z']), (tree_sat3['ssfr'] * 1e11), '-r')
axins.set_xlim(8.5, 12)
axins.set_ylim(-0.02 * 10 ** (-11) * 1e11, 2 * 10 ** (-11) * 1e11)
axins.vlines(t_quench3, min(tree_sat3['ssfr']), max(tree_sat3['ssfr'] + 0.5), colors='black', linestyles='--', lw=1,
             zorder=10)
axins.hlines(10 ** (-10.93) * 1e11, tmin, tmax, colors='black', linestyles='--', lw=1, zorder=10)
axins.tick_params(axis='x', direction='in', length=6, width=1.2, which='major', colors='black',
                  grid_color='black', grid_alpha=0.4)
axins.tick_params(axis='y', which='major', direction='in', length=6, width=0.9, colors='black',
                  grid_color='black', grid_alpha=0.4)

axins.text(10.95, +1.45 * 10 ** (-11) * 1e11, r'$ \mathtt{SSFRcut}  $',
           rotation=0,
           horizontalalignment='center',
           verticalalignment='top',
           multialignment='center',
           fontdict=small_font)
axins.set_yticks([0, 1])
axins.set_yticklabels(('0', r'$\mathtt{1\times 10^{-11}}$'))
axins.set_xlim(8, 12)
ax4.set_xlabel('Time since Big Bang [Gyr]')
ax3.set_ylabel(r'$\mathtt{SSFR \ \ [yr^{-1}]}$')
plt.subplots_adjust(hspace=0, wspace=0.33)

# Do the plot code
# plt.savefig('ssfr_dist.pdf', format='pdf', dpi=1200, pad_inches=0.1, bbox_inches='tight')
plt.show()

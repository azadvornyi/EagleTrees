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


t= time.time()
sim100 = 'RefL0100N1504'
sim25 = 'RefL0025N0376'

sim_box = int(argv[8])
#============
# sim = int(argv[3])  # set simulation here
#============

if sim_box == 100:
    sim = sim100
    box_size = 100
    host_mass = 1E13
else:
    sim = sim25
    box_size = 25
    host_mass = 1E12

con = sql.connect('twh176', password='tt408GS2')






gid = int(argv[1])
fof = int(argv[2])
sub = int(argv[3])
dx = float(argv[4])
dy = float(argv[5])
dz = float(argv[6])
sat_index = int(argv[7])  # select satellite here
host_index = gid
#print(sim, gid,dx,dy,dz,box_size,sat_index)
# satellite selection around a host at z = 0. Hosts close to the wall of the box are handled as well
normal_sat_query= 'SELECT \
             S.GalaxyID as Sgid, \
             H.GalaxyID as Hgid, \
             FOF.Group_R_Crit200 as r_vir \
             FROM \
             {0}_SubHalo as H, \
             {0}_SubHalo as S, \
             {0}_FOF as FOF \
             WHERE \
             H.GalaxyID = {1:.0f} \
             and H.GroupID = FOF.GroupID \
             and 0.0033*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_x  - (S.CentreOfPotential_x - FLOOR((S.CentreOfPotential_x+{2:.0f})/ {5:.0f}))) \
             and 0.0033*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_y  - (S.CentreOfPotential_y - FLOOR((S.CentreOfPotential_y+{3:.0f})/ {5:.0f}))) \
             and 0.0033*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_z  - (S.CentreOfPotential_z - FLOOR((S.CentreOfPotential_z+{4:.0f})/ {5:.0f}))) \
             and S.Snapnum = 28 \
             and S.MassType_Star between 1E9 and 1E12'.format(sim, gid, dx, dy, dz, box_size)

sats_info = sql.execute_query(con, normal_sat_query)

sat_len = len(sats_info["Sgid"])
#print(sat_len)
# retrieving SubGroupNumber and FoF for a satellite
fof_sub_query = \
    'SELECT \
S.GalaxyID as gid, \
S.GroupNumber as fof, \
S.SubGroupNumber as sub \
FROM \
{0}_SubHalo as S \
WHERE \
S.GalaxyID = {1}'.format(sim, sats_info["Sgid"][sat_index])

fof_sub_info = sql.execute_query(con, fof_sub_query)

# print(fof_sub_info,'fof subinfo')


# Host tree query
tree_host_query = 'SELECT \
             DES.GalaxyID as gid, \
             DES.TopLeafID as tlid, \
             PROG.Redshift as z, \
             PROG.MassType_DM as mdm, \
             PROG.MassType_Star as ms,\
             AP.SFR / (AP.Mass_Star+0.0001) as ssfr, \
             PROG.CentreOfPotential_x as copx, \
             PROG.CentreOfPotential_y as copy, \
             PROG.CentreOfPotential_z as copz \
         FROM \
             {0}_Subhalo as PROG, \
             {0}_Subhalo as DES, \
             {0}_Aperture as AP \
         WHERE \
             DES.SnapNum = 28 \
             and DES.GroupNumber = {1:.0f} \
             and DES.SubGroupNumber = {2:.0f} \
             and PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID \
             and AP.ApertureSize = 100 \
             and AP.GalaxyID = PROG.GalaxyID \
         ORDER BY \
             DES.MassType_Star desc, \
             PROG.Redshift asc, \
             PROG.MassType_Star desc'.format(sim, fof, sub)



# print(sats_info['Sgid'][1])
# print(int(fof_sub_info['fof']),'fof')

#test_tree = sql.execute_query(con,tree_query_generator(sim1,fof_sub_test['fof'], fof_sub_test['sub'] ))
#
#print(test_tree)

# Satellite tree query
tree_sat_query = 'SELECT \
             DES.GalaxyID as gid, \
             DES.TopLeafID as tlid,\
             PROG.Redshift as z, \
             PROG.MassType_DM as mdm, \
             PROG.MassType_Star as ms, \
             AP.SFR / (AP.Mass_Star+0.0001) as ssfr, \
             PROG.CentreOfPotential_x as copx, \
             PROG.CentreOfPotential_y as copy, \
             PROG.CentreOfPotential_z as copz \
         FROM \
             {0}_Subhalo as PROG, \
             {0}_Subhalo as DES, \
             {0}_Aperture as AP \
         WHERE \
             DES.SnapNum = 28 \
             and DES.GroupNumber = {1:.0f} \
             and DES.SubGroupNumber = {2:.0f} \
             and PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID \
             and AP.ApertureSize = 100 \
             and AP.GalaxyID = PROG.GalaxyID \
         ORDER BY \
             DES.MassType_Star desc, \
             PROG.Redshift asc, \
             PROG.MassType_Star desc'.format(sim, int(fof_sub_info['fof']), int(fof_sub_info['sub']))


#building host and satellite trees


try:
    with open("sim{0}/{1}/host_{1}".format(sim100,host_index), 'rb') as f3:
        tree_host = pickle.load(f3)
    #tree_host = np.load("sim{0}/{1}/host_{1}.npy".format(sim25,host_index))
except FileNotFoundError:
    tree_host = sql.execute_query(con, tree_host_query)


tree_sat = sql.execute_query(con, tree_sat_query)

data = tree_sat

#print(np.sqrt((tree_host['copx'][0] - tree_sat['copx'][0])**2 + (tree_host['copy'][0] - tree_sat['copy'][0])**2 + (tree_host['copz'][0] - tree_sat['copz'][0])**2))
#print((tree_host['copx'][0] - tree_sat['copx'][0]))


# retrieving virial mass of the host
host_r_vir_query = 'SELECT \
             PROG.Redshift as z,\
             FOF.Group_R_Crit200 as r_vir\
         FROM \
             {0}_Subhalo as PROG, \
             {0}_Subhalo as DES, \
             {0}_Aperture as AP,\
             {0}_FOF as FOF \
         WHERE \
             DES.SnapNum = 28 \
             and FOF.GroupID = PROG.GroupID \
             and DES.GroupNumber = {1:.0f} \
             and DES.SubGroupNumber = {2:.0f} \
             and PROG.GalaxyID between DES.GalaxyID and DES.TopLeafID \
             and AP.ApertureSize = 100 \
             and AP.GalaxyID = PROG.GalaxyID \
         ORDER BY \
             DES.MassType_Star desc, \
             PROG.Redshift asc, \
             PROG.MassType_Star desc'.format(sim, fof, sub)


#host_r_vir = sql.execute_query(con, host_r_vir_query)

try:
    with open("sim{0}/{1}/host_r_vir_{1}".format(sim100,host_index), 'rb') as f5:
        host_r_vir = pickle.load(f5)
    #tree_host = np.load("sim{0}/{1}/host_{1}.npy".format(sim25,host_index))
except FileNotFoundError:
    host_r_vir = sql.execute_query(con, host_r_vir_query)




# converting z to Gyr
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
        if distances[counter] < 0.0033 * r_vir[counter]:
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

# print(t_infall,"infall time")
# print(t_quench,'t_quench')
# print(tiq,'t_quench -  t_infall, Gyr')
# # print(data['ms'][0],'M_solar @ z = 0')

sat_mass = data['ms'][0]

# checking if the host became a satellite of itself
if (tree_host['copx'][0] - tree_sat['copx'][0]) == 0:
    pass
elif (t_quench or t_infall) == 0 :
    pass
else:
    f = open("data_plot_{0}.txt".format(sim100), "a")
    f.write("{0:.0f} {1:.0f} {2} {3} {4} {5} {6}\n".format(host_index, sat_index, tree_host['mdm'][0], sat_mass,
                                                 t_infall, t_quench, data['ssfr'][0]))
    f.close()


this_sat_id = sat_index

dirName = 'sim{1}/{0}'.format(host_index,sim100)

# Create target directory & all intermediate directories if don't exists

if (tree_host['copx'][0] - tree_sat['copx'][0]) == 0:
    pass
# elif (t_quench or t_infall) == 0 :
#     pass
else:
    try:
        os.makedirs(dirName)
        print("Directory ", dirName, " Created ")
        with open("sim{1}/{0}/host_{0}".format(host_index, sim100), 'wb') as f0:
            pickle.dump(tree_host, f0, pickle.HIGHEST_PROTOCOL)
        with open("sim{1}/{0}/host_r_vir_{0}".format(host_index, sim100), 'wb') as f1:
            pickle.dump(host_r_vir, f1, pickle.HIGHEST_PROTOCOL)
    except FileExistsError:
        pass


    with open("sim{2}/{0}/sat_{1}".format(host_index,sat_index, sim100), 'wb') as f2:
        pickle.dump(tree_sat, f2, pickle.HIGHEST_PROTOCOL)




print(time.time() - t,"sec ")

# plt.plot(time_z, radius)
# plt.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0025,c = 'red')
# plt.ylabel('radius, (Mpc)')
# plt.hlines(first_approach_r,0,max(time_z))
# plt.vlines(t_infall,0, 1)
# plt.savefig('dist.pdf', format ='pdf')
# plt.clf()
#
# # #print(data['sns'])
# plt.plot(data['z'], data['copx'], '-r')
# plt.plot(data['z'], data['copy'], '-g')
# plt.plot(data['z'], data['copz'], '-b')
# plt.savefig('cop.pdf')
#
# plt.clf()
# plt.semilogy(times_Gyr(data['z']), data['ms'], '-k')
# plt.savefig('ms.pdf', format='pdf')
#
# plt.clf()
# plt.semilogy(times_Gyr(data['z']), data['ssfr']*1e9, '-r')
# plt.vlines(t_quench,0, 1)
# plt.savefig('ssfr.pdf', format='pdf')

#plt.show()
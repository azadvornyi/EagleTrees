import eagleSqlTools as sql
import numpy as np
from sys import argv
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import time
from scipy import interpolate
import pickle
from sys import argv
import sys
import os
sim100 = 'RefL0100N1504'
sim25 = 'RefL0025N0376'
fof1 = 2
sub1 = 7
#============
sim = sim25  # set simulation here
#============

if sim == sim100:
    box_size = 100
    host_mass = 1E13
else:
    box_size = 25
    host_mass = 1E12

con = sql.connect('twh176', password='tt408GS2')

host_index = 2
#int(argv[1])  # select host here
sat_index = 1 #int(argv[2])   # select satellite here



#  selecting hosts above given mass
host_query = "SELECT \
             SH.GalaxyID as gid, \
             SH.GroupNumber as fof,\
             SH.SubGroupNumber as sub, \
             SH.CentreOfPotential_x as copx, \
             SH.CentreOfPotential_y as copy, \
             SH.CentreOfPotential_z as copz, \
             FOF.Group_R_Crit200 as r_vir \
             FROM\
             {0}_SubHalo as SH, \
             {0}_FOF as FOF\
             WHERE \
             FOF.GroupID = SH.GroupID \
             and SH.SnapNum = 28 \
             and SH.MassType_DM > {1}".format(sim, host_mass)

host_ids = sql.execute_query(con, host_query)

print(len(host_ids))

if int(host_ids['sub'][host_index]) != 0:
    sys.exit()

#print((host_ids['gid'][0]))
host_id = 13892596
host_id1 = 152519
t = time.time()

he_index = host_index

t1 = time.time()
# Handling hosts that are close to the sides of the box
def edge_proximity(host_info,box_size): # copy this
    offsetx = np.array([])
    offsety = np.array([])
    offsetz = np.array([])
    gid = np.array([])
    fof = np.array([])
    sub = np.array([])
    halfbox = box_size / 2
    for host in host_info:
        x, y, z = 0, 0, 0
        r_vir = host['r_vir'] * 0.0025
        if (host['copx'] - r_vir < 0):
            x = halfbox
        elif (host['copx'] + r_vir > box_size):
            x = -halfbox
        elif (host['copy'] - r_vir < 0):
            y = halfbox
        elif (host['copy'] + r_vir > box_size):
            x = -halfbox
        elif (host['copz'] - r_vir < 0):
            z = halfbox
        elif (host['copz'] + r_vir > box_size):
            z = -halfbox

        gid = np.append(gid, host['gid'])
        fof = np.append(fof, host['fof'])
        sub = np.append(sub, host['sub'])
        offsetx = np.append(offsetx, x)
        offsety = np.append(offsety, y)
        offsetz = np.append(offsetz, z)

    edge_hosts = np.core.records.fromarrays([gid, fof, sub,offsetx, offsety, offsetz], names='gid,fof,sub,dx,dy,dz',
                                                formats='i4, i4, i4,float64,float64,float64')

    return edge_hosts

host_on_edge = edge_proximity(host_ids, box_size)# copy this

print(host_on_edge,"this is host on edge")
# print(host_ids)
print(time.time() - t1, "s, moving to edge")
t1 = time.time()
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
             and 0.0025*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_x  - (S.CentreOfPotential_x - FLOOR((S.CentreOfPotential_x+ {2:.0f})/ {5:.0f}))) \
             and 0.0025*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_y  - (S.CentreOfPotential_y - FLOOR((S.CentreOfPotential_y+ {3:.0f})/ {5:.0f}))) \
             and 0.0025*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_z  - (S.CentreOfPotential_z - FLOOR((S.CentreOfPotential_z+ {4:.0f})/ {5:.0f}))) \
             and S.Snapnum = 28 \
             and S.MassType_Star between 1E9 and 1E12'.format(sim, host_on_edge['gid'][he_index], host_on_edge['dx'][he_index],
                                                              host_on_edge['dy'][he_index],host_on_edge['dz'][he_index],box_size)
sats_info = sql.execute_query(con, normal_sat_query)
#print(sats_info)
print( len(sats_info["Sgid"]))
print(time.time() - t1, "s, sat selection around a host")
t1 =time.time()
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

print(fof_sub_info,'fof subinfo of a astellite')
print(time.time() - t1, "s, retrieving SubGroupNumber and FoF for a satellite")

t1 = time.time()
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
             PROG.MassType_Star desc'.format(sim, int(host_ids['fof'][host_index]), int(host_ids['sub'][host_index]))



try:
    with open('sim25/{0}/host_{0}'.format(host_index), 'rb') as f3:
        tree_host = pickle.load(f3)
    #tree_host = np.load("sim{0}/{1}/host_{1}.npy".format(sim25,host_index))
except FileNotFoundError:
    print("calculating the host tree")
    tree_host = sql.execute_query(con, tree_host_query)



print(time.time() - t1, "s, Host tree query")
# print(sats_info['Sgid'][1])
# print(int(fof_sub_info['fof']),'fof')

#test_tree = sql.execute_query(con,tree_query_generator(sim1,fof_sub_test['fof'], fof_sub_test['sub'] ))
#
#print(test_tree)
t1 = time.time()
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

tree_sat = sql.execute_query(con, tree_sat_query)
print(time.time() - t1, "s, sat tree query")
data = tree_sat


#print(np.sqrt((tree_host['copx'][0] - tree_sat['copx'][0])**2 + (tree_host['copy'][0] - tree_sat['copy'][0])**2 + (tree_host['copz'][0] - tree_sat['copz'][0])**2))
#print((tree_host['copx'][0] - tree_sat['copx'][0]))

t1 = time.time()
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
             PROG.MassType_Star desc'.format(sim, int(host_ids['fof'][host_index]), int(host_ids['sub'][host_index]))



# host_r_vir = sql.execute_query(con, host_r_vir_query)
# print(len(host_r_vir))
try:
    with open("sim25/{0}/host_r_vir_{0}".format(host_index), 'rb') as f4:
        host_r_vir = pickle.load(f4)
except FileNotFoundError:
    print("calculating R_vir")
    host_r_vir = sql.execute_query(con, host_r_vir_query)




print(time.time() - t1, "s, R-vir host")
print("Host coordinates")
print(tree_host['copx'][0],tree_host['copy'][0],tree_host['copz'][0])
print("="*20)
print("Sat coordinates")
print(tree_sat['copx'][0],tree_sat['copy'][0],tree_sat['copz'][0])
print("R_vir")
print(sats_info['r_vir'][0])
print(host_r_vir['r_vir'])
print(len(host_r_vir['r_vir']))
print(len(tree_sat))
#print(host_ids['r_vir'][17])

#print(sats_info)

t1 = time.time()
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
        if distances[counter] < 0.0025 * r_vir[counter]:
            first_approach = np.append(first_approach,distances[counter])
            t_infall = np.append(t_infall, time_[counter])
    distances = list(reversed(distances))
    time_ = list(reversed(time_))
    time_ = times_Gyr(time_)

    return distances, time_


#print(tree_sat['ssfr'],'this is sat')

radius, time_z = moving_to_origin(tree_host, tree_sat, box_size, host_r_vir['r_vir'])

print(time.time() - t, 'it took this many seconds')


# calculating t_quench - t_infall
def t_infall_and_quench(r_vir, dist, time_dist, ssfr, gen_time, sat_time):
    r_vir_func = interpolate.interp1d(gen_time, r_vir )
    dist_func = interpolate.interp1d(time_dist, dist )

    precise_time_d = np.linspace(min(time_dist),max(time_dist), 150)
    precise_time = np.linspace(min(gen_time), max(gen_time), 150)

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
    interp_ssfr = ssfr_func(precise_time)
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
                                                                times_Gyr(tree_host['z']), times_Gyr((tree_sat['z'])))

print(time.time() - t1, "s, calculating distances")

# print(t_infall,"infall time")
# print(t_quench,'t_quench')
print(tiq,'t_quench -  t_infall, Gyr')
print(data['ms'][0],'M_solar @ z = 0')

dirName = 'sim25/{0}'.format(host_index)

# Create target directory & all intermediate directories if don't exists
try:
    os.makedirs(dirName)
    print("Directory ", dirName, " Created ")
    with open("sim25/{0}/host_{0}".format(host_index), 'wb') as f0:
        pickle.dump(tree_host, f0, pickle.HIGHEST_PROTOCOL)

except FileExistsError:
    print("Directory ", dirName, " already exists")


with open("sim25/{0}/sat_{1}".format(host_index,sat_index), 'wb') as f2:
    pickle.dump(tree_sat, f2, pickle.HIGHEST_PROTOCOL)

with open("sim25/{0}/host_r_vir_{0}".format(host_index), 'wb') as f1:
    pickle.dump(host_r_vir, f1, pickle.HIGHEST_PROTOCOL)

print(t_infall,'t_infall')
print(t_quench,"t_quench")


plt.plot(time_z, radius)
plt.plot(times_Gyr(host_r_vir['z']), host_r_vir['r_vir']*0.0025,c = 'red')
plt.ylabel('radius, (Mpc)')
plt.hlines(first_approach_r,0,max(time_z))
plt.vlines(t_infall,0, 1)
plt.savefig('dist1.pdf', format ='pdf')
plt.clf()

# #print(data['sns']) tree_host
plt.plot(data['z'], data['copx'], '-r')
plt.plot(data['z'], data['copy'], '-g')
plt.plot(data['z'], data['copz'], '-b')
plt.plot(tree_host['z'], tree_host['copx'], '-o')
plt.plot(tree_host['z'], tree_host['copy'], 'olive')
plt.plot(tree_host['z'], tree_host['copz'], 'cyan')
plt.savefig('cop.pdf')

plt.clf()
plt.semilogy(times_Gyr(data['z']), data['ms'], '-k')
plt.savefig('ms.pdf', format='pdf')

plt.clf()
plt.semilogy(times_Gyr(data['z']), data['ssfr']*1e9, '-r')
plt.vlines(t_quench,0, 1)
plt.savefig('ssfr.pdf', format='pdf')

#plt.show()
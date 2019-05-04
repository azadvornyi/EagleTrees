
from sys import argv

import eagleSqlTools as sql
import numpy as np
from sys import argv
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import time

from sys import argv

sim100 = 'RefL0100N1504'
sim25 = 'RefL0025N0376'

#============
sim_box = int(argv[2])
#============
sim = int(argv[2])

if sim_box == 100:
    sim = sim100
    box_size = 100
    host_mass = 1E13
else:
    sim = sim25
    box_size = 25
    host_mass = 1E12

con = sql.connect('twh176', password='tt408GS2')

host_index = int(argv[1])  # argv[1] select host here
#print(argv[1])

sat_index = 0   # select satellite here



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
#print(len(host_ids))
#print((host_ids['gid'][0])
t = time.time()
he_index = host_index

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
             and 0.0025*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_x  - (S.CentreOfPotential_x - FLOOR((S.CentreOfPotential_x+{2:.0f})/ {5:.0f}))) \
             and 0.0025*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_y  - (S.CentreOfPotential_y - FLOOR((S.CentreOfPotential_y+{3:.0f})/ {5:.0f}))) \
             and 0.0025*FOF.Group_R_Crit200 > ABS( H.CentreOfPotential_z  - (S.CentreOfPotential_z - FLOOR((S.CentreOfPotential_z+{4:.0f})/ {5:.0f}))) \
             and S.Snapnum = 28 \
             and S.MassType_Star between 1E9 and 1E12'.format(sim, host_on_edge['gid'][he_index], host_on_edge['dx'][he_index],
                                                              host_on_edge['dy'][he_index],host_on_edge['dz'][he_index],box_size)
sats_info = sql.execute_query(con, normal_sat_query)

#print( len(sats_info["Sgid"]))


try:
    sat_len = len(sats_info["Sgid"])
    f = open("hostid_satnumber_{0}.txt".format(sim_box), "a")
    f.write("{0:.0f} {1:.0f} \n".format(host_index, sat_len ))
    f.close()
    #print(sat_len)
except TypeError as error:
    f = open("hostid_satnumber_{0}.txt".format(sim_box), "a")
    f.write("{0:.0f} 0 \n".format(host_index))
    f.close()






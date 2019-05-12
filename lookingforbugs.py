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
sim25 = 'RecalL0025N0752'

#============
sim_box = int(argv[1])
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
        r_vir = host['r_vir'] * 0.0033
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

host_on_edge = edge_proximity(host_ids, box_size)

for i in range(len(host_ids["gid"])):
    if host_ids["sub"][i] == 0:
        try:
            f = open("hostid_rvir_33_{0}.txt".format(sim_box), "a")
            f.write("{0} {1} {2} {3} {4} {5} \n".format(host_on_edge["gid"][i],host_on_edge["fof"][i],host_on_edge["sub"][i],
                                                        host_on_edge["dx"][i],host_on_edge["dy"][i],host_on_edge["dz"][i]))
            f.close()
            #print(sat_len)
        except TypeError as error:
            pass









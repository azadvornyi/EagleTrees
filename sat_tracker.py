
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
sim_box = int(argv[7])
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

host_index = int(argv[7])  # argv[1] select host here
#print(argv[1])

sat_index = 0   # select satellite here



gid = int(argv[1])
fof = int(argv[2])
sub = int(argv[3])
dx = float(argv[4])
dy = float(argv[5])
dz = float(argv[6])





# Handling hosts that are close to the sides of the box
# copy this


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
             and S.MassType_Star between 1E9 and 1E12'.format(sim, gid, dx, dy, dz, box_size)


#print( len(sats_info["Sgid"]))


sats_info = sql.execute_query(con, normal_sat_query)
print(sats_info)
# try:
#     f = open("hostid_satnumber_{0}.txt".format(sim_box), "a")
#     f.write("{0} {1} {2} {3} {4} {5} {6}\n".format(gid, fof, sub , dx, dy, dz, len(sats_info) ))
#     f.close()
#     print(len(sats_info))
# except TypeError as error:
#     pass



# f = open("hostid_satnumber_{0}.txt".format(sim_box), "a")
    # f.write("{0:.0f} 0 \n".format(host_index))
    # f.close()






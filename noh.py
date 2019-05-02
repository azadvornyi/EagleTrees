
from sys import argv

import eagleSqlTools as sql
import numpy as np
from sys import argv
import matplotlib as mpl
mpl.use('TkAgg')

sim100 = 'RefL0100N1504'
sim25 = 'RefL0025N0376'

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
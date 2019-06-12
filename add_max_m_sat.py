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

host_id = "19054212"


sim = "RefL0100N1504"
sat_id = 7



box_size = 100


host_gid,sat_num,M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504_with_v.txt", unpack=True, usecols=(0,1,2, 3, 4, 5, 6))
# M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RecalL0100N0752.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))



keep_these = np.where(t_q != 14.134423953091941)  #
and_these = np.where((t_i <= t_q))
q_before_i = np.where((t_i >= t_q))
unquenched = np.where(t_q == 14.134423953091941)

# print(m_s[keep_these])
a = set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
c = np.array(list(c))





for i in c:

    with open("sim{0}_v/{1}/host_r_vir_{1}".format(sim, host_id), 'rb') as f1:
        host_r_vir = pickle.load(f1)

    with open("sim{0}_v/{1}/host_{1}".format(sim, int(host_gid[i])), 'rb') as f2:
        tree_host = pickle.load(f2)

    with open("sim{0}_v/{1}/sat_{2}".format(sim, int(host_gid[i]), int(sat_num[i])), 'rb') as f3:
        tree_sat = pickle.load(f3)
    with open("sim{0}_v/{1}/sat_vel_{2}".format(sim, int(host_gid[i]), int(sat_num[i])), 'rb') as f4:
        tree_sat_v = pickle.load(f4)


    f = open("Data_plot_100_with_max_mass".format(sim), "a")
    f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(int(host_gid[i]),int(sat_num[i]),M_h[i],
                                                                    max(tree_sat["ms"]), t_i[i], t_q[i],
                                                                    t_q[i]-t_i[i], sfr[i]))
    f.close()






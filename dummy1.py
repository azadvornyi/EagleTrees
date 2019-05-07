# from sys import argv
#
# a = int(argv[1])
#
# print(a)
#print('-')
# import pickle
# #
# # with open("sim25/0/host_r_vir_0", 'rb') as f4:
# #     host_r_vir = pickle.load(f4)
# #
# # print(len(host_r_vir['r_vir']))
import numpy as np
n = np.loadtxt("hostid_satnumber_100.txt", unpack=True, usecols=(6,))
print(np.sum(n))
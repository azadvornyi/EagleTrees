# from sys import argv
#
# a = int(argv[1])
#
# print(a)
#print('-')
#with open('trunk.pickle', 'rb') as f:
#    t = pickle.load(f)
import numpy as np

a = np.load("sim25/0/sat_0.npy")

print(a["copy"])
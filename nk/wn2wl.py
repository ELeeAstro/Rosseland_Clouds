import numpy as np


data = np.loadtxt('NH3_ori.dat',delimiter=',')
wn = data[:,0]
n = data[:,1]
k = data[:,2]

nwl = len(wn)

wl = 1.0/wn * 1e4

wl = wl[::-1]
n = n[::-1]
k = k[::-1]

for l in range(nwl):
    print(wl[l],n[l],k[l])

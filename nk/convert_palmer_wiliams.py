import numpy as np
import matplotlib.pylab as plt

fname = 'palmer_williams_wl.dat'
data = np.loadtxt(fname)


wl = data[:,0]
pw_75_n = data[:,1]
pw_75_k = data[:,2]
pw_85_n = data[:,3]
pw_85_k = data[:,4]
pw_96_n = data[:,5]
pw_96_k = data[:,6]

for i in range(len(wl)):
  print(wl[i],pw_96_n[i],pw_96_k[i])


exit()

fig = plt.figure()

plt.plot(wl[:-1],pw_96_n[:-1])
plt.plot(wl[:-1],pw_96_k[:-1])


plt.xscale('log')

plt.show()

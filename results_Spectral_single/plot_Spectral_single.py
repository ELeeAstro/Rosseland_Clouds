import numpy as np
import matplotlib.pylab as plt

sp = 'MgSiO3_amorph'

# read in k_ext data
fname = sp+'_single.txt'

f = open(fname, 'r')
line = f.readline().split()
r_med = float(line[0])
nwl = int(line[1])

data = np.loadtxt(fname,skiprows=1)
wl = data[:,0]
kext = data[:,1]
a = data[:,2]
g = data[:,3]

fig = plt.figure()

plt.plot(wl,kext)
plt.xscale('log')
plt.yscale('log')

plt.ylabel(r'k$_{\rm ext}$ [cm$^{2}$]')
plt.xlabel(r'wavelength [$\mu$m]')

plt.title('r: ' + str(r_med))

plt.savefig(sp+'_kext.png',dpi=144,bbox_inches='tight')


fig = plt.figure()

plt.plot(wl,a)
plt.xscale('log')

plt.title('r: ' + str(r_med) )

plt.ylabel(r'$\omega$ [-]')
plt.xlabel(r'wavelength [$\mu$m]')

plt.savefig(sp+'_a.png',dpi=144,bbox_inches='tight')

fig = plt.figure()

plt.plot(wl,g)
plt.xscale('log')

plt.title('r: ' + str(r_med))

plt.ylabel(r'$g$')
plt.xlabel(r'wavelength [$\mu$m]')

plt.savefig(sp+'_g.png',dpi=144,bbox_inches='tight')


plt.show()
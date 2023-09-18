import numpy as np
import matplotlib.pylab as plt

## Species names to plot ##
sp = 'Mg2SiO4_amorph'
## ----- ##

dirc = 'results_lognorm/'

fname = dirc+sp+'_ln.txt'
f = open(fname,'r')
head = f.readline().split()
r_med = float(head[0])
sigma = float(head[1])
rmin = float(head[2])
rmax = float(head[3])
iint = int(head[4])
nwl = int(head[5])

print(r_med, sigma, rmin, rmax, iint, nwl)

data = np.loadtxt(fname,skiprows=1)

wl = data[:,0]
k_norm = data[:,1]
a = data[:,2]
g = data[:,3]

# Calculate normalized size distribution for plotting
rad = np.logspace(np.log10(rmin),np.log10(rmax),iint)
nd_dist = np.zeros(iint)
for r in range(iint):
  nd_dist[r]= (1.0/(rad[r] * np.sqrt(2.0*np.pi) * np.log(sigma))) \
  * np.exp(-(np.log(rad[r]/r_med))**2/(2.0 * np.log(sigma)**2))


fig = plt.figure()

# Calculate some log-normal distribution properties here
r_mean = r_med * np.exp(0.5 * np.log(sigma)**2)
r_mode = np.exp(np.log(r_med) - np.log(sigma)**2)
r_eff = r_med * np.exp(5.0/2.0 * np.log(sigma)**2)

plt.title(str(r_med) + ' ' + str(sigma))

plt.plot(rad,nd_dist,c='black',lw=2)

plt.axvline(r_med,c='darkgreen',label='r$_{med}$')
plt.axvline(r_mean,c='darkred',label='r$_{mean}$')
plt.axvline(r_mode,c='darkmagenta',label='r$_{mode}$')
plt.axvline(r_eff,c='darkcyan',label='r$_{eff}$')

plt.xscale('log')
plt.yscale('log')

plt.ylabel('n$_{d, norm}$ [um$^{-1}$]')
plt.xlabel('radius [$\mu$m]')

plt.legend()

plt.savefig('plots_lognorm/'+sp+'_dist_norm.png',dpi=144,bbox_inches='tight')

fig = plt.figure()

plt.title(str(r_med) + ' ' + str(sigma))

plt.plot(wl,k_norm,c='black',lw=2)

plt.xscale('log')
plt.yscale('log')

plt.ylabel('k$_{norm}$ [cm$^{2}$ particle$^{-1}$]')
plt.xlabel('wavelength [$\mu$m]')

plt.savefig('plots_lognorm/'+sp+'_k_norm.png',dpi=144,bbox_inches='tight')

fig = plt.figure()

plt.title(str(r_med) + ' ' + str(sigma))

plt.plot(wl,a,c='black',lw=2)

plt.xscale('log')

plt.ylabel('a')
plt.xlabel('wavelength [$\mu$m]')

plt.savefig('plots_lognorm/'+sp+'_a.png',dpi=144,bbox_inches='tight')

fig = plt.figure()

plt.title(str(r_med) + ' ' + str(sigma))

plt.plot(wl,g,c='black',lw=2)

plt.xscale('log')

plt.ylabel('g')
plt.xlabel('wavelength [$\mu$m]')

plt.savefig('plots_lognorm/'+sp+'_g.png',dpi=144,bbox_inches='tight')

plt.show()

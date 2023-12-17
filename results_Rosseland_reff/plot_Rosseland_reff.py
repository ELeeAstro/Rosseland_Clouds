import numpy as np
import matplotlib.pylab as plt

sp = 'MgSiO3_amorph'

# read in k_ext data
fname = sp+'_kext.txt'

f = open(fname, 'r')
line = f.readline().split()
na = int(line[0])
nT = int(line[1])
sigma = float(line[2])
line = f.readline().split()
a = np.zeros(na)
a[:] = line[:]
line = f.readline().split()
reff = np.zeros(na)
reff[:] = line[:]
line = f.readline().split()
T= np.zeros(nT)
T[:] = line[:]

print(na, nT)
print(a[:])
print(reff[:])
print(T[:])

kext = np.log10(np.loadtxt(fname,skiprows=4))


fig = plt.figure()

kext_lev = np.linspace(np.min(kext),np.max(kext),10)

con_kext = plt.contourf(T,reff,kext,kext_lev)
plt.colorbar(con_kext,label=r'$\log_{10}$ k$_{\rm ext}$ [cm$^{2}$]')
plt.yscale('log')

plt.ylabel('effective radius [um]')
plt.xlabel('temperature [K]')


plt.title('sigma: ' + str(sigma))

plt.savefig(sp+'_kext.png',dpi=144,bbox_inches='tight')

# read in a data
fname = sp+'_a.txt'

f = open(fname, 'r')
line = f.readline().split()
na = int(line[0])
nT = int(line[1])
sigma = float(line[2])
line = f.readline().split()
a = np.zeros(na)
a[:] = line[:]
line = f.readline().split()
reff = np.zeros(na)
reff[:] = line[:]
line = f.readline().split()
T= np.zeros(nT)
T[:] = line[:]

print(na, nT)
print(a[:])
print(reff[:])
print(T[:])

w = np.loadtxt(fname,skiprows=4)

fig = plt.figure()

w_lev = np.linspace(np.min(w),np.max(w),10)

con_w = plt.contourf(T,reff,w,w_lev)
plt.colorbar(con_w,label=r'$\omega$ [-]')
plt.yscale('log')

plt.ylabel('effective radius [um]')
plt.xlabel('temperature [K]')

plt.title('sigma: ' +  str(sigma))

plt.savefig(sp+'_a.png',dpi=144,bbox_inches='tight')

# read in g data
fname = sp+'_g.txt'

f = open(fname, 'r')
line = f.readline().split()
na = int(line[0])
nT = int(line[1])
sigma = float(line[2])
line = f.readline().split()
a = np.zeros(na)
a[:] = line[:]
line = f.readline().split()
reff = np.zeros(na)
reff[:] = line[:]
line = f.readline().split()
T= np.zeros(nT)
T[:] = line[:]

print(na, nT)
print(a[:])
print(reff[:])
print(T[:])

g = np.loadtxt(fname,skiprows=4)

fig = plt.figure()

g_lev = np.linspace(np.min(g),np.max(g),10)

con_g = plt.contourf(T,reff,g,g_lev)
plt.colorbar(con_g,label=r'$g$ [-]')
plt.yscale('log')

plt.ylabel('effective radius [um]')
plt.xlabel('temperature [K]')

plt.title('sigma: ' +  str(sigma))

plt.savefig(sp+'_g.png',dpi=144,bbox_inches='tight')


plt.show()
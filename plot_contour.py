import numpy as np
import matplotlib.pylab as plt


## Species names to plot ##
sp = 'MgSiO3'
## ----- ##

nt = 22
nr = 51

fname = 'rosselandMean_RTtable.txt'
f = open(fname,'r')
dim = f.readline().split()
rad = np.array(f.readline().split())
rad = rad.astype(np.float) 
temp = np.array(f.readline().split())
temp = temp.astype(np.float)

data = np.loadtxt('rosselandMean_RTtable.txt',skiprows=1)
rad = data[0,:]
temp = data[1,:]

print(rad)
print(temp)

Qext = np.loadtxt('results/'+sp+'_rosselandMean_qext.txt',skiprows=1)
Qsca = np.loadtxt('results/'+sp+'_rosselandMean_qscat.txt',skiprows=1)
gg = np.loadtxt('results/'+sp+'_rosselandMean_gg.txt',skiprows=1)

Qext_lev = np.linspace(0,np.max(Qext),10)
Qsca_lev = np.linspace(0,np.max(Qsca),10)
gg_lev = np.linspace(np.min(gg),np.max(gg),10)


fig = plt.figure()

con_Qext = plt.contourf(temp,rad,Qext,Qext_lev)
plt.colorbar(con_Qext,label='Q$_{ext}$')
plt.yscale('log')

plt.ylabel('radius [m]')
plt.xlabel('temperature [K]')

plt.savefig('plots/'+sp+'_Qext.png',dpi=144,bbox_inches='tight')


fig = plt.figure()

con_Qsca = plt.contourf(temp,rad,Qsca,Qsca_lev)
plt.colorbar(con_Qsca,label='Q$_{sca}$')
plt.yscale('log')

plt.ylabel('radius [m]')
plt.xlabel('temperature [K]')

plt.savefig('plots/'+sp+'_Qsca.png',dpi=144,bbox_inches='tight')

fig = plt.figure()

con_gg = plt.contourf(temp,rad,gg,gg_lev)
plt.colorbar(con_gg,label='g')
plt.yscale('log')

plt.ylabel('radius [m]')
plt.xlabel('temperature [K]')

plt.savefig('plots/'+sp+'_g.png',dpi=144,bbox_inches='tight')

plt.show()

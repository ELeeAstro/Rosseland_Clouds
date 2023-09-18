import numpy as np
import matplotlib.pylab as plt


## Species names to plot ##
sp = 'Mg2SiO4_amorph'
## ----- ##

fname = 'RTtable.txt'
f = open(fname,'r')
dim = f.readline().split()
rad = np.array(f.readline().split())
rad = rad.astype(float)
temp = np.array(f.readline().split())
temp = temp.astype(float)

print(dim)
print(len(rad),rad)
print(len(temp),temp)

Qext = np.loadtxt('results_Ross/'+sp+'_rosselandMean_qext.txt',skiprows=1)
Qsca = np.loadtxt('results_Ross/'+sp+'_rosselandMean_qscat.txt',skiprows=1)
gg = np.loadtxt('results_Ross/'+sp+'_rosselandMean_gg.txt',skiprows=1)

print(Qext.shape)

Qext_lev = np.linspace(0,np.max(Qext),10)
Qsca_lev = np.linspace(0,np.max(Qsca),10)
gg_lev = np.linspace(np.min(gg),np.max(gg),10)


fig = plt.figure()

con_Qext = plt.contourf(temp,rad,Qext,Qext_lev)
plt.colorbar(con_Qext,label='Q$_{ext}$')
plt.yscale('log')

plt.ylabel('radius [um]')
plt.xlabel('temperature [K]')

plt.savefig('plots_Ross/'+sp+'_Qext.png',dpi=144,bbox_inches='tight')


fig = plt.figure()

con_Qsca = plt.contourf(temp,rad,Qsca,Qsca_lev)
plt.colorbar(con_Qsca,label='Q$_{sca}$')
plt.yscale('log')

plt.ylabel('radius [um]')
plt.xlabel('temperature [K]')

plt.savefig('plots_Ross/'+sp+'_Qsca.png',dpi=144,bbox_inches='tight')

fig = plt.figure()

con_gg = plt.contourf(temp,rad,gg,gg_lev)
plt.colorbar(con_gg,label='g')
plt.yscale('log')

plt.ylabel('radius [um]')
plt.xlabel('temperature [K]')

plt.savefig('plots_Ross/'+sp+'_g.png',dpi=144,bbox_inches='tight')

plt.show()

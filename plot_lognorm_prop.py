import numpy as np
import matplotlib.pylab as plt

# Calculate some log-normal distribution properties given a median value and sigma value

r_med = 1.0
sig = 2.0

r_mode = r_med * np.exp(-np.log(sig)**2)
r_mean = r_med * np.exp(1.0/2.0 * np.log(sig)**2)
r_eff = r_med * np.exp(5.0/2.0 * np.log(sig)**2)
r_v = r_med * np.exp(7.0/2.0 * np.log(sig)**2)

# Ackerman & Marley (2001) - characteristic velocity radius
fsed = 2.0
alpha = 1.3
r_w = r_med * fsed**(alpha) * np.exp((alpha + 6.0)/2.0 * np.log(sig)**2)

print('sigma: ', sig)
print('r_med: ', r_med)
print('r_mode: ', r_mode)
print('r_mean: ', r_mean)
print('r_eff: ', r_eff)
print('r_v: ', r_v)

print('---')
print('fsed: ', fsed)
print('alpha: ', alpha)
print('r_w: ', r_w)

# Calculate normalized (or given total number density) number density distribution

# Size distribution sizes
na = 100
amin = 1e-3
amax = 1e2
a = np.logspace(np.log10(amin),np.log10(amax),na) * 1e-4 # Convert to cm

N0 = 1.0 # N0 = 1 is the normalised version
nd = np.zeros(na)

for i in range(na):
  nd[i] = N0/(a[i]*np.sqrt(2.0*np.pi)*np.log(sig)) * np.exp(-(np.log(a[i]/(r_med*1e-4))**2)/(2.0*np.log(sig)**2))

# Calculate specific nd for calculated a values
nd_mode = N0/(r_mode*1e-4*np.sqrt(2.0*np.pi)*np.log(sig)) * np.exp(-(np.log(r_mode/r_med)**2)/(2.0*np.log(sig)**2))
nd_mean = N0/(r_mean*1e-4*np.sqrt(2.0*np.pi)*np.log(sig)) * np.exp(-(np.log(r_mean/r_med)**2)/(2.0*np.log(sig)**2))
nd_eff = N0/(r_eff*1e-4*np.sqrt(2.0*np.pi)*np.log(sig)) * np.exp(-(np.log(r_eff/r_med)**2)/(2.0*np.log(sig)**2))
nd_v = N0/(r_v*1e-4*np.sqrt(2.0*np.pi)*np.log(sig)) * np.exp(-(np.log(r_v/r_med)**2)/(2.0*np.log(sig)**2))
nd_w = N0/(r_w*1e-4*np.sqrt(2.0*np.pi)*np.log(sig)) * np.exp(-(np.log(r_w/r_med)**2)/(2.0*np.log(sig)**2))

a[:] = a[:] * 1e4 # Convert back to um

fig = plt.figure

plt.plot(a,nd,lw=2,c='grey')

plt.vlines(r_mode,0,nd_mode,colors='darkred',lw=2,label='Mode',ls='dashed')
plt.vlines(r_mean,0,nd_mean,colors='darkgreen',lw=2,label='Mean',ls='dashed')
plt.vlines(r_eff,0,nd_eff,colors='darkmagenta',lw=2,label='Eff.',ls='dashed')
plt.vlines(r_v,0,nd_v,colors='darkcyan',lw=2,label='Vol.',ls='dashed')
plt.vlines(r_w,0,nd_w,colors='goldenrod',lw=2,label='Char.',ls='dashed')

plt.xlabel('Particle size [$\mu$m]')
plt.ylabel('Size distribution [cm$^{-3}$ cm$^{-1}$]')

plt.legend()

plt.xscale('log')
plt.yscale('log')

plt.show()
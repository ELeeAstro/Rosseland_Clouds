import numpy as np

na = 51
amin = 1e-8
amax = 1e-4

Tmin = 100.0
Tmax = 2500.0
dT = 50.0 

a = np.logspace(np.log10(amin),np.log10(amax),na)
T = np.arange(Tmin,Tmax+dT,dT)

#print(len(a),a)

#print(len(T), T)


print(len(a), len(T))
print(" ".join(str(l) for l in a[:]))
print(" ".join(str(l) for l in T[:]))

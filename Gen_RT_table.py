import numpy as np

# Number of radius points
na = 51

# min and max radius [um]
amin = 1e-1
amax = 100.0

# min and max temperature grid
Tmin = 500.0
Tmax = 2500.0

# Temperature spacing 
dT = 50.0 

a = np.logspace(np.log10(amin),np.log10(amax),na)
T = np.arange(Tmin,Tmax+dT,dT)

#print(len(a),a)

#print(len(T), T)

# Print to terminal, can use '>' pipe to send output to a file of choice
print(len(a), len(T))
print(" ".join(str(l) for l in a[:]))
print(" ".join(str(l) for l in T[:]))

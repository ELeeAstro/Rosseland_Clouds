import numpy as np

# Number of wavelengths
nwl = 10000

# min and max wavelengths [um]
wlmin = 0.2
wlmax = 1000.0

# Set up wavelength grid (logspaced recommended)
wl = np.logspace(np.log10(wlmin),np.log10(wlmax),nwl)

# For GCM spectral models at certain bands, apply the wavelengths manually
# e.g. below we have the 11 bands of Tiffany Kataria of MITgcm fame.
#nwl = 12
#wl = [0.260, 0.420, 0.610, 0.850, 1.320, 2.020, 2.500, 3.500, 4.400, 8.70, 20.00 ,324.68]

# Output wavelengths to file
output_name = 'wavelengths.txt'
output = open(output_name,'w')
output.write(str(nwl) + '\n')
for i in range(nwl):
  output.write(str(wl[i]) + '\n')

## Now we do the radius and temperature spacing

# Number of radius points
na = 51

# min and max radius [um]
amin = 1e-2
amax = 20.0

# size grid (logspaced recommended)
a = np.logspace(np.log10(amin),np.log10(amax),na)

# Number of temperature points
nT = 21

# min and max temperature grid
Tmin = 200.0
Tmax = 2000.0

T = np.linspace(Tmin,Tmax,nT)

# Output radius and T grid to file
output_name = 'RTtable.txt'
output = open(output_name,'w')
output.write(str(na) + ' ' +  str(nT) + '\n')
output.write(" ".join(str(l) for l in a[:]) + '\n')
output.write(" ".join(str(l) for l in T[:]) + '\n')

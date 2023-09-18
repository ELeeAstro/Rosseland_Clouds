import numpy as np

# Number of wavelengths
nwl = 10000

# min and max wavelengths [um]
wlmin = 0.2
wlmax = 1000.0

# Set up wavelength grid
wl = np.logspace(np.log10(wlmin),np.log10(wlmax),nwl)

# Output wavelengths to file
output_name = 'wavelengths_RM.txt'
output = open(output_name,'w')
output.write(str(nwl) + '\n')
for i in range(nwl):
  output.write(str(wl[i]) + '\n')

import numpy as np

nwl = 10000
wl = np.logspace(np.log10(0.1),np.log10(1000),nwl)


output_name = 'wavelengths.txt'
output = open(output_name,'w')
output.write(str(nwl) + '\n')
for i in range(nwl):
  output.write(str(wl[i]) + '\n')

# Rosseland_Clouds
 Rosseland mean calculator for exoplanet cloud species

 Warning: This code is under active development and undergoing a refresh of an older code for easier use. The README is not complete!

 ## How to use

 The Gen_input.py file will generate the input files wavelengths.txt and RTtable.txt - the wavelengths, radii and temperatures of the calculations are contained here.
 You can set things like the wavelength grid to perform the calculations on and temperature range of the Rosseland mean calculations.

 The Rosseland_clouds.nml file contains input for the fortran code to read in. Here you set things like the path to the input files, species and nk data to use. For log-normal calculations the distribution parameters can be set here.

 To compile the code, enter 'make' in the src directory. This will compile the code and make the executable 'Rosseland_clouds' in the main directory. Then run the code by running ./Rosseland_clouds in the main directory.

 ## Modes: 

 ### 1. Rosseland mean with single particle size

 This takes the particle size array in RTtable as a single particle size (delta function), and calculates the Rosseland mean normalised opacity, single scattering albedo and asymmetry factor.

 ### 2. Rosseland mean with lognormal size distribution

 ### 3. Rosseland mean with lognormal reff size (similar to single sizes)

  This takes the particle size array in RTtable as the median particle size of a log-normal distribution, then calculates the effective particle size given the sigma value in the namelist.
  Then it calculates the Rosseland mean normalised opacity, single scattering albedo and asymmetry factor similar to the single particle size case.

 ### 4. Spectral with single particle sizes at lognormal effective size

 This takes the particle size array and wavelength array and calculates the effective particle size of the log-normal parameters.
 It then calculates the normalised kext, a and g at those specific wavelengths at the effective particle size. 
 This mode is useful for building tables for spectral RT models in the GCM.

 ## Output

 Output is generated into specific results_* directories for your method and species. We in

 ## Example use in GCMs

 One question to ask is, how do I actually use this data inside the GCM or RT code I have?

 ## Useful stuff

 The plot_log_norm_prop.py python code can calculate and plot log-normal distribution properties given a median and sigma value (edit this in the python script).
 This is useful to visualise how the size-distribution looks before calculation, as well as a useful calculator for useful values such as the effective radius and mass weighted radius etc.

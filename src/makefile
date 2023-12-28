# for gfortran Compiler
#===================

F90 = gfortran
F90LINKER =  gfortran

#Debugging and development flags
#FFLAGS	= -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
#FFLAGS	= -Og -pipe -Wall -Wextra -g -fbacktrace

#Serial flags
FFLAGS  = -O3 -pipe

#Parallel flags
#FFLAGS	= -O3 -pipe -fopenmp

# for ifort Compiler
#====================

#F90 = ifort
#F90LINKER = ifort

#FFLAGS   = -O0 -g -traceback -fpp -prec-div -fp-model source -fpe0 -ipo
#FFLAGS    = -O0 -g -traceback -xHost -fpp -fp-model source -qopenmp -ipo
#FFLAGS    = -O3 -xHost -qopenmp -fpp -fp-model source -ipo

#====================

DEFS      =
INCLUDES  =
LFLAGS    = $(FFLAGS)


OBJECTS = \
lxmie_mod.o \
Rosseland_data_mod.o \
Rosseland_single.o \
Rosseland_lognorm.o \
Rosseland_reff.o \
Spectral_reff.o \
Rosseland_main.o

# executable statement
EXECS  = ../Rosseland_clouds

.SUFFIXES : .o .f90 .f

default: Rosseland_clouds

Rosseland_clouds: $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

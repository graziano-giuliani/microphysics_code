
FC = gfortran
#FCFLAGS = -I. `nf-config --fflags` -O3
FCFLAGS = -I. `nf-config --fflags` -ffpe-trap=invalid,zero,overflow \
          -O0 -g -fcheck=all -fbacktrace
LIBS = `nf-config --flibs`

.SUFFIXES: .F90 .f90 .o

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

all :: microphysics_code

OBJS = mod_constants.o mod_dimensions.o mod_micro_common.o mod_micro_nogtom.o

microphysics_code: microphysics_code.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(OBJS) $(LIBS)

mod_constants.o: mod_constants.F90
mod_dimensions.o: mod_dimensions.F90 mod_constants.o
mod_micro_common.o: mod_micro_common.F90 mod_constants.o mod_dimensions.o
mod_micro_nogtom.o: mod_micro_nogtom.F90 mod_constants.o mod_dimensions.o mod_micro_common.o

clean:
	rm -f *.mod *.o microphysics_code

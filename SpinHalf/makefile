#FC = gfortran
#FFLAGS = -O3 -Wall -Wextra -g -fbacktrace -fcheck=all -pedantic -cpp #-pg
#-ffpe-trap= invalid, zero, overflow
#-Werror -ffree-line-length-500 #-fpp -D IFORT -qopenmp 
#-pg -> per debugging   $ gprof <program-name> <gmon.out>
#PFLAGS = -fopenmp 
#LIBS = -llapack -lblas

FC = ifort
FFLAGS = -O2 -fpp -heap-arrays -D IFORT -g -warn all -traceback -check bounds -pedantic #,noexternal, -check all
PFLAGS = -qopenmp 
LIBS = -mkl -Wl,-rpath=/opt/pgi-8/common/intel/lib/intel64_lin/ 

SRC = mod_print.f90 mod_functions.f90 mod_obs.f90 mod_states.f90 mod_exp.f90 mod_mat.f90 sorts.f90 mod_entanglement.f90 #mataid.f90 expokit.f90 exp_sparse.f90
MOD = ${SRC:.f90=.o} #substitute .f90 with .o in all 'SRC' files

#This is the first command, which is always executed. If any '.f90' files in 'SRC' have been modified, it compiles them. 
%.o : %.f90 #creation of all *.o files DEPENDS on *.f90
	$(FC) $(FFLAGS) $(PFLAGS) -c $<

dummy: $(MOD) dummy.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)
dummy2: $(MOD) dummy2.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

spectrum_swap_LR: $(MOD) spectrum_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

solvable_spectrum_swap_LR: $(MOD) solvable_spectrum_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

spectrum_Khemani: $(MOD) spectrum_Khemani.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

correlations_swap_LR: $(MOD) correlations_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

dynamics_swap_LR: $(MOD) dynamics_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

decay_swap_LR: $(MOD) decay_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

decay_kick0_swap_LR: $(MOD) decay_kick0_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

comparison_kick0_skipping_swap_LR: $(MOD) comparison_kick0_skipping_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

all: spectrum_swap_LR, correlations_swap_LR, dynamics_swap_LR, decay_swap_LR


clean:
	@mv -f *.o *.mod ../Trash

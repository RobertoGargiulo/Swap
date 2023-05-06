FC =  gfortran
FFLAGS = -O3 -Wall -Wextra -g -fbacktrace -fcheck=all -pedantic #-pg
#-ffpe-trap= invalid, zero, overflow
#-Werror -ffree-line-length-500 #-fpp -D IFORT -qopenmp 
#-pg -> per debugging   $ gprof <program-name> <gmon.out>
SRC = mod_print.f90 mod_functions.f90 mod_obs.f90 mod_mat.f90 mod_states.f90 mod_exp.f90 sorts.f90 mod_entanglement.f90 #mataid.f90 expokit.f90 exp_sparse.f90
LIBS = -llapack -lblas
PFLAGS = -fopenmp #-qopenmp
MOD = ${SRC:.f90=.o} #substitute .f90 with .o

#This is the first command, which is always executed. If any '.f90' files have been modified, it compiles them. 
%.o : %.f90 #creation of all *.o files DEPENDS on *.f90
	$(FC) $(FFLAGS) $(PFLAGS) -c $<

swap: $(MOD) swap.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

swap_dense_Sz0: $(MOD) swap_dense_Sz0.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

swap_spectrum_Sz0: $(MOD) swap_spectrum_Sz0.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

swap_decay: $(MOD) swap_decay.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

swap_entanglement_Sz0: $(MOD) swap_entanglement_Sz0.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

swap_corr_z_Sz0: $(MOD) swap_corr_z_Sz0.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

swap_overlap: $(MOD) swap_overlap.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

swap_degeneracies_Sz0: $(MOD) swap_degeneracies_Sz0.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

spectrum_swap_LR: $(MOD) spectrum_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

solvable_spectrum_swap_LR: $(MOD) solvable_spectrum_swap_LR.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

spectrum_flip: $(MOD) spectrum_flip.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

all: swap swap_dense_Sz0 swap_spectrum_Sz0 swap_decay swap_entanglement_Sz0 #dense_Sz0 sparse_Sz0 spectrum_Sz0


clean:
	@mv -f *.o *.mod ../../Trash

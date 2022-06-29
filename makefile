FC = gfortran
FFLAGS = -O3 -Wall -Wextra -fbacktrace -fcheck=all #-fsyntax-only 
SRC = mod_print.f90 mod_genmat.f90 mod_exp.f90 mataid.f90 expokit.f90 exp_sparse.f90 
LIBS = -llapack -lblas
PFLAGS = -fopenmp
MOD = ${SRC:.f90=.o} #substitute .f90 with .o

#This is the first command, which is always executed. If any '.f90' files have been modified, it compiles them. 
%.o : %.f90 #creation of all *.o files DEPENDS on *.f90
	$(FC) $(FFLAGS) $(PFLAGS) -c $<

swap: $(MOD) swap.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

prova: $(MOD) prova.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

prova2: $(MOD) prova2.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

prova3: $(MOD) prova3.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

mag_avg: mag_avg.f90
	$(FC) $(FFLAGS) -o $@ $@.f90

all: swap mag_avg prova prova2 prova3


clean:
	@mv -f *.o *.mod ../Trash

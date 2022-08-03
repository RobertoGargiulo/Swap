FC =  gfortran
FFLAGS = -O3 -Wall -Wextra -g -fbacktrace -fcheck=all -fopenmp -pedantic 
#-ffpe-trap= invalid, zero, overflow
#-Werror -ffree-line-length-500 #-fpp -D IFORT -qopenmp 
#-pg -> per debugging   $ gprof <program-name> <gmon.out>
SRC = mod_print.f90 mod_genmat.f90 mod_exp.f90 mataid.f90 expokit.f90 exp_sparse.f90 
LIBS = -llapack -lblas
PFLAGS = #-qopenmp #-fopenmp
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

prova4: $(MOD) prova4.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

prova5: $(MOD) prova5.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

dense_Sz0: $(MOD) dense_Sz0.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

sparse_Sz0: $(MOD) sparse_Sz0.o
	$(FC) $(FFLAGS) $(PFLAGS) -o $@ $(MOD) $@.o $(LIBS)

mag_avg: mag_avg.f90
	$(FC) $(FFLAGS) -o $@ $@.f90

min_time_avg: $(MOD) min_time_avg.o
	$(FC) $(FFLAGS) -o $@ $(MOD) $@.o $(LIBS)

all: swap mag_avg prova prova2 prova3 min_time_avg


clean:
	@mv -f *.o *.mod ../Trash

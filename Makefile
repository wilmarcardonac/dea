FC     = gfortran
LC     = $(FC)
EXE    = dea
INCDIR = /usr/local/include/fgsl
LIBDIR = /usr/local/lib
F_FL   = -O3 -Wall -I$(INCDIR) -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -g
LIB_FL = -L$(LIBDIR) -lfgsl 

#####################
OBJ   = fiducial.o input.o background.o initial_conditions.o perturbations.o diff_asde.o radau5.o decsol.o dc_decsol.o 

def:	$(OBJ) $(OBJNR) $(OBJODE)
	$(LC) $(F_FL) $(OBJ) $(OBJNR) $(OBJODE) -o $(EXE)  $(LIB_FL)

%.o:	%.f90
	$(FC) $(F_FL) -c $<

%.o:	%.F90
	$(FC) $(F_FL) -c $<

%.o:	%.f
	$(FC) $(F_FL) -c $<

clean :
	rm -f *.o *.mod *.stb *~ *.il $(EXE)

### put dependencies here ###

diff_asde.o :	diff_asde.f90 fiducial.o input.o background.o initial_conditions.o perturbations.o   

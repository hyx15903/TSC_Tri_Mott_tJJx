.SUFFIX:
.SUFFIXES: .o .mod .F .f .F90 .f90
F90FLAGS= -g -O2 -xhost -qopenmp -qmkl=parallel -traceback -mcmodel=medium 
FC= ifort 
DIR= $1
TAG= tjmodel
OBJS=common.o block.o measure.o  get_sys.o den_mat.o  eigen.o dmrg_sub.o main.o
TAG = tjmodel
all : 	$(TAG)
$(TAG):	$(OBJS) 
	$(FC)  $(OBJS) -qopenmp -qmkl=parallel -o $(TAG)  
.f90.o:
	$(FC) -c $< $(F90FLAGS)
clean: 
	rm -f *.o *.mod core 

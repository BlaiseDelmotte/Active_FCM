FC	= ifort
LD	= $(FC)
F90	= $(FC)
SWP	= 
RM	= /bin/rm -f
MP	=
#OLEVEL	= -O3 -mnoieee-fp -unroll -msse2
#OLEVEL	= -O3 -msse2 -unroll
OLEVEL	= -O3
FOPTS	= -r8 
#FOPTS	= -r8 -ftrapuv -check bounds
FFLAGS	= $(FOPTS) $(OLEVEL)
FFTWDIR = /usr/local/P3DFFT
LIBS    = $(FFTWDIR)/lib/libp3dfft.a $(FFTWDIR)/lib/libfftw3.a -lmpi
INCLUDE =-I$(FFTWDIR)/include

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=



MODOBJS	= \
module_var.o \
MPI_structures.o


FLUOBJS	= \
get_dims.o \
utilities_cpu.o \
meshing.o \
interp_lag0.o \
interp_lag1.o \
interp_lag2.o \
interp_lag3.o \
interp_lag4.o \
main.o


PAROBJS	= 

SCLOBJS	= 


MOD	= *.mod
OBJS	= $(MODOBJS) $(FLUOBJS) $(PAROBJS) $(SCLOBJS) 
EXEC    = /users/p1208/pfede/SPLITFIELD/EXEC/splitfield

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) $(INCLUDE) -o $@

clean:
	$(RM) $(OBJS) $(MOD) $(EXEC)

cpart:
	$(RM) $(PAROBJS) $(MOD) $(EXEC)

.SUFFIXES: .o .f90 

.f90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<

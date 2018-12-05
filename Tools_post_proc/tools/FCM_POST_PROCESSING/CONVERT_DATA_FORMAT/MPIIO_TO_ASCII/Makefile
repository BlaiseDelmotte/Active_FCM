FC	= mpif90
LD	= $(FC)
F90	= $(FC)
SWP	= 
RM	= /bin/rm -f
MP	=
#OLEVEL	= -O3 -mnoieee-fp -unroll -msse2
OLEVEL	= -O3 -unroll
FOPTS	= -r8 -lstdc++
#FOPTS	= -r8 -ftrapuv -check bounds -g
FFLAGS	= $(FOPTS) $(OLEVEL)

FFTWDIR = /PRODCOM/FFTW/fftw-3.3.1-intel64v2011.7.256/glibc-2.11
P3DFFTDIR = /PRODCOM/FFTW/p3dfft.2.5.1-intel64v2011.7.256/glibc-2.11
TECPLOTDIR = /PRODCOM/TECPLOT/tec360-2011r2/x86-64-2.9
LIBS	= $(P3DFFTDIR)/lib/libp3dfft.a $(FFTWDIR)/lib/libfftw3.a $(TECPLOTDIR)/lib/tecio.a
INCLUDE = -I$(P3DFFTDIR)/include -I$(FFTWDIR)/include -I$(TECPLOTDIR)/include 

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=


MODOBJS	= \



FLUOBJS	= \
read_MPIIO.o\
read_var_MPIIO.o\
mpiio_to_ascii_pos_only.o\
#mpiio_to_ascii_pos_noper_only.o\
#mpiio_to_ascii.o\


FCMOBJS	= \



PAROBJS	= \


SCLOBJS	= \



MOD	= *.mod
OBJS	= $(MODOBJS) $(FLUOBJS) $(FCMOBJS) $(PAROBJS) $(SCLOBJS) 
EXEC    = mpi2ascii

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) $(INCLUDE) -o $@

clean:
	$(RM) $(OBJS) $(MOD) $(EXEC)

cpart:
	$(RM) $(PAROBJS) $(MOD) $(EXEC)

.SUFFIXES: .o .f90 

.f90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<

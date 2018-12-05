FC	= mpif90
LD	= $(FC)
F90	= $(FC)
SWP	= 
RM	= /bin/rm -f
MP	=
#OLEVEL	= -O3 -mnoieee-fp -unroll -msse2
OLEVEL	= -O3 -unroll
FOPTS	= -r8 
#FOPTS	= -r8 -ftrapuv -check bounds -g
FFLAGS	= $(FOPTS) $(OLEVEL)

FFTWDIR = /PRODCOM/FFTW/fftw-3.3.1-intel64v2011.7.256/glibc-2.11
P3DFFTDIR = /PRODCOM/FFTW/p3dfft.2.5.1-intel64v2011.7.256/glibc-2.11
LIBS	= $(P3DFFTDIR)/lib/libp3dfft.a $(FFTWDIR)/lib/libfftw3.a 
INCLUDE = -I$(P3DFFTDIR)/include -I$(FFTWDIR)/include 

LDFLAGS	= $(FOPTS) $(OLEVEL)
PROF	=


MODOBJS	= \



FLUOBJS	= \
postreat_fluid.o\
read_MPIIO.o

FCMOBJS	= \



PAROBJS	= \


SCLOBJS	= \



MOD	= *.mod
OBJS	= $(MODOBJS) $(FLUOBJS) $(FCMOBJS) $(PAROBJS) $(SCLOBJS) 
EXEC    = postreat_fluid

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) $(INCLUDE) -o $@

clean:
	$(RM) $(OBJS) $(MOD) $(EXEC)

cpart:
	$(RM) $(PAROBJS) $(MOD) $(EXEC)

.SUFFIXES: .o .f90 

.f90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<

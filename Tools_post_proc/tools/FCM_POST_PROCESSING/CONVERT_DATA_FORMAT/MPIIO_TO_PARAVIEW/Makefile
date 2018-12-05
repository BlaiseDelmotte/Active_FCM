include Makefile.inc


MODOBJS	= \



FLUOBJS	= \
ell_mat_vectors.o\
print_paraview_part_ell.o\
read_MPIIO.o\
read_var_MPIIO.o\
mpiio_to_paraview.o\
print_paraview_part.o\
#~ mpiio_to_paraview_pos_only.o\
#~ print_paraview_part_pos_only.o\


FCMOBJS	= \



PAROBJS	= \


SCLOBJS	= \



MOD	= *.mod
OBJS	= $(MODOBJS) $(FLUOBJS) $(FCMOBJS) $(PAROBJS) $(SCLOBJS) 
EXEC    = mpi2para

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) $(INCLUDE) -o $@

clean:
	$(RM) $(OBJS) $(MOD) $(EXEC)

cpart:
	$(RM) $(PAROBJS) $(MOD) $(EXEC)

.SUFFIXES: .o .f90 

.f90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<

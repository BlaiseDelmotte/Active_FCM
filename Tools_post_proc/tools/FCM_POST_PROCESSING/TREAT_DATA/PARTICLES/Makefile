include Makefile.inc


MODOBJS	= \
mod_bucket_variable.o\



FLUOBJS	= \
bucket_allocation.o\
bucket_distribution.o\
bucket_build_map.o\
bucket_inumber.o\
compute_bucket_lists.o\
compute_connectivity.o\
fit_lin_log.o\
height_distrib.o\
mean_vel.o\
msd_disp.o\
msd_only.o\
nematic_order.o\
pdf_velocities.o\
pdf_rotations.o\
pdf_orientations.o\
polar_order.o\
polar_order_only.o\
position_velocity_orient_correl.o\
position_velocity_orient_correl_bidisperse.o\
position_velocity_orient_correl_prolate.o\
rdf_only.o\
read_MPIIO.o\
read_var_MPIIO.o\
rot_mat_time.o\
rot_mat_time_vectors.o\
short_time_self_diff.o\
vel_autocorrel.o\
treat_particles_bidisperse.o\
#treat_particles.o\
#treat_pos_only_percolation.o\
# treat_msd_rdf_only.o\
#treat_pos_only.o\
#treat_particles.o\
#treat_fluctu_only.o\








FCMOBJS	= \



PAROBJS	= \


SCLOBJS	= \



MOD	= *.mod
OBJS	= $(MODOBJS) $(FLUOBJS) $(FCMOBJS) $(PAROBJS) $(SCLOBJS) 
EXEC    = treat_part

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) $(INCLUDE) -o $@

clean:
	$(RM) $(OBJS) $(MOD) $(EXEC)

cpart:
	$(RM) $(PAROBJS) $(MOD) $(EXEC)

.SUFFIXES: .o .f90 

.f90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<

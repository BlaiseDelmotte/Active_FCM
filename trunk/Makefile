include Makefile.inc


MODOBJS	= \
mod_dns_dim.o \
mod_check_cpu.o \
mod_cputime.o \
mod_ensight_var.o \
mod_fcm_bucket_variable.o \
mod_fcm_forcing_variable.o \
mod_fcm_part_variable.o \
mod_fluid_variable.o \
mod_forcing.o \
mod_geometric_variable.o \
mod_param_phys.o\
mod_particle_parallel.o \
mod_rhs_variable.o \
mod_scalar_variable.o \
mod_statistics.o \
mod_work_arrays.o \
mod_collision.o \
mod_random_var.o \
MPI_structures.o \
utilities_cpu.o \
exchange_part.o


FLUOBJS	= \
add_forcing2.o \
adv_fluid.o\
allocate_arrays.o\
build_rhs_vort.o \
check_cfl.o \
cputime_info.o \
corrtwopts.o \
divergence.o \
fluid_gradient.o \
fluid_prediction.o \
inforun.o\
init_run.o \
initiation_fluid.o\
initiation_forcing.o\
meshing.o \
openclose.o \
print_fluid.o \
proj_divfree.o \
print_timestat.o \
print_lastcycle.o \
readparam.o \
save_fluid.o\
save_forcing.o\
spec3d.o\
stat_fluid.o \
tremain.o \
vorticity.o \
main.o

FCMOBJS	= \
fcm_add_monopole_forcing.o \
fcm_add_rotlet_forcing.o \
fcm_add_stresslet_forcing.o \
fcm_add_swimming_forcing.o \
fcm_adv_particle_position.o \
fcm_adv_particle_orientation.o \
fcm_adv_particle_orientation_pswim_only.o \
fcm_algorithm.o \
fcm_algorithm_brownian.o \
fcm_barrier_force.o \
fcm_barrier_force_no_bucket.o \
fcm_bucket_distribution.o \
fcm_bucket_allocations.o \
fcm_bucket_build_map.o \
fcm_bucket_inumber.o \
fcm_compute_bucket_lists.o \
fcm_compute_gaussian_env.o \
fcm_compute_gaussian_support.o \
fcm_compute_gaussian_env_function.o \
fcm_compute_gaussian_dip_env_function.o \
fcm_compute_gaussian_grad_dip_env_function.o \
fcm_compute_gaussian_grad_env_function.o \
fcm_compute_gaussian_coeff_potdip_env_function.o \
fcm_compute_gaussian_env_sq.o \
fcm_compute_self_ros_tensor.o \
fcm_compute_self_vel_tensor.o \
fcm_conjugate_gradient_stresslet.o \
fcm_distrib_monopole.o \
fcm_distrib_rotlet.o \
fcm_distrib_squirming_forcing.o \
fcm_distrib_stresslet.o \
fcm_fluid_prediction.o \
fcm_fluid_prediction_brownian.o \
fcm_fluid_rate_of_strain.o \
fcm_gb_shifted_uniaxial.o \
fcm_initiation.o\
fcm_initiation_barrier_magnitude.o \
fcm_initiation_particle_gaussian.o \
fcm_initiation_particle_orientation.o \
fcm_initiation_particle_position.o \
fcm_initiation_particle_pos_orient.o \
fcm_initiation_particle_velocity.o \
fcm_initiation_phase_shift.o \
fcm_initiation_rep_wall_magnitude.o \
fcm_inv_3_3.o \
fcm_kinematic_contact_handling_no_bucket.o\
fcm_kinematic_contact_handling_bucket_mono.o\
fcm_kinematic_contact_handling_bucket_mono_arg.o\
fcm_KS_second_step.o\
fcm_l2norm.o \
fcm_mirror_forces_x.o \
fcm_oseen_dipole.o \
fcm_oseen_quadrupole.o \
fcm_prod_scal.o \
fcm_random_gaussian_variable.o \
fcm_rate_of_strain_filter.o \
fcm_readparam.o \
fcm_remove_self_ros.o \
fcm_remove_self_vel.o \
fcm_repulsive_wall.o \
fcm_RF_setup.o\
fcm_rfd_second_step.o\
fcm_rotation_filter.o \
fcm_RS_setup.o\
fcm_RS_setup_symmetrize.o\
fcm_save_particle_characteristics.o \
fcm_save_particle_kinematics.o \
fcm_save_particle_dynamics.o \
fcm_sum_forces_RS_fourier.o \
fcm_time_stroke.o \
fcm_velocity_filter.o \
fcm_write_info.o \
fcm_zero_field_variable.o\
fcm_zero_part_variable.o\
#fcm_kinematic_contact_handling_para_arg.o\


PAROBJS	= \
adv_particle_position.o \
adv_particle_velocity.o \
boundary_particle.o\
collision_detection.o\
collision_interaction.o\
dgqt.o\
destsv.o \
ell_full2eig.o\
ell_full2low.o\
ell_low2chol.o\
ell_pair_closest_points.o\
ell_pair_intersec.o\
ell_pair_separate.o\
ellu_pt_near_far.o\
ell_line_proj.o\
filtering_fluid.o\
gaspart_transfer.o \
initiation_particle_position.o \
initiation_particle_velocity.o \
initiation_particle_scalar.o \
init_ellipsoid.o \
init_ellipsoid_2D.o \
interpH.o \
interp_lag1.o\
interp_lag2.o\
interp_lag3.o\
interp_lag4.o \
intersec_ellipsoid.o \
collision.o \
particle_tracking.o \
print_lagfunction.o\
seed_ellipsoid.o \
seed_ellipsoid_2D.o \
stat_particle.o \
save_particle.o \
save_particle_time.o

SCLOBJS	= \
initiation_scalar.o\
adv_scalar.o\
stat_scalar.o\
save_scalar.o\
spec3d_scalar.o


MOD	= *.mod
OBJS	= $(MODOBJS) $(FLUOBJS) $(FCMOBJS) $(PAROBJS) $(SCLOBJS) 

EXEC    = ns3d


$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) $(LIBS) $(INCLUDE) -o $@

clean:
	$(RM) $(OBJS) $(MOD) $(EXEC)

cpart:
	$(RM) $(PAROBJS) $(MOD) $(EXEC)

.SUFFIXES: .o .f90 .F90 

.f90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<

.F90.o:
	$(FC)  -c $(FFLAGS) $(INCLUDE) $<

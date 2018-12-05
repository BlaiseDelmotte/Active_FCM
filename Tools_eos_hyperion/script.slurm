#!/bin/bash
#SBATCH --nodes=13
#SBATCH --ntasks=260
#SBATCH --ntasks-per-node=20
#SBATCH --threads-per-core=1
#SBATCH -J NS3D_SW
#SBATCH --time=150:00:00
#SBATCH --exclusive
#SBATCH --mail-user=bdelmott@imft.fr



export WALLTIME=540000 # walltime pour PBS en secondes


###################################
# CPUFREQ #
###################################
export FREQ=2800000
#Available(KHz): 0(2801000Turbo) 1(2800000) 2(2600000) 3(2400000) 4(2200000)
# 5(2000000) 6(1800000) 7(1600000) 8(1400000) 9(1200000)
clush -bw $SLURM_NODELIST "sudo /home_nfs/script/admin/pwrb f $FREQ"

###################################
# COMPILER #
###################################
# IntelMPI
module purge
module load intel/14.0.2.144
module load intelmpi/4.1.3.049
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

###################################
# Profiling MPI #
###################################
#module load mpiprof_sp1.1.106_mpi_4.1.0.030 # pour intelmpi
###################################
#

###################################
NNODES=$SLURM_NNODES
cores=$SLURM_NTASKS
###################################
# OMP_NUM_THREADS #
###################################
export OMP_NUM_THREADS=1
echo DIR_$FREQ\_$cores\_$OMP_NUM_THREADS\_$SLURM_JOBID
#export EXE=$SLURM_SUBMIT_DIR/ns3d
export EXE=$HOME/NS3D_active_brownian_bidisperse/trunk_run_tumble/ns3d
cd $SLURM_SUBMIT_DIR
echo '------------------------------------------------'
echo ' NODE USED = '$SLURM_NODELIST
echo ' SLURM_JOBID = '$SLURM_JOBID
echo ' OMP_NUM_THREADS = '$OMP_NUM_THREADS
echo ' CORES = '$cores
echo '------------------------------------------------'
echo 'ldd '$EXE
ldd $EXE
export mpi_tasks=$cores
export NODE_TYPE=HTOFF
time srun --cpu_bind=map_cpu:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 -n $cores $EXE $WALLTIME >  ns3d_$cores\_$SLURM_JOBID.log

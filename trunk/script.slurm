#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH -J NS3D_SW
#SBATCH -C HSW24
#SBATCH --time=01:00:00
#SBATCH --mail-user=blaise.delmotte@ladhyx.polytechnique.fr

export WALLTIME=3600 # walltime pour PBS en secondes


###################################
# COMPILER USED FOR FFW3 and P3DFFT #
###################################
# OpenMPI
module purge
module load intel/17.2 openmpi/intel/2.0.2

###################################
# Profiling MPI #
###################################
#module load mpiprof_sp1.1.106_mpi_4.1.0.030 # pour intelmpi
###################################

###################################
# OUTPUT_DIRECTORY #
###################################
export dirname=$SCRATCHDIR/$SLURM_JOBID
#export dirname=$HOME/$SLURM_JOBID
echo '------------------------------------------------'
echo ' JOB WILL BE RUN FROM : ' $dirname
echo '------------------------------------------------'


export postprocdir=$HOME/NS3D_active_brownian_bidisperse/Tools_post_proc/tools/FCM_POST_PROCESSING/TREAT_DATA/PARTICLES/
echo '------------------------------------------------'
echo ' PARTICLE POST_PROC EXECUTABLE COPIED FROM : ' $postprocdir
echo '------------------------------------------------'


mkdir $dirname
cp ns3d $dirname
cp *.in $dirname
cp *.ini $dirname
cp $postprocdir/treat_part $dirname
cd $dirname
ls -al

export outname=$STOREDIR/results_$SLURM_JOBID
mkdir $outname
echo '------------------------------------------------'
echo ' DATA WILL BE SAVED IN : ' $outname
echo '------------------------------------------------'

###################################
NNODES=$SLURM_NNODES
cores=$SLURM_NTASKS
###################################
# OMP_NUM_THREADS #
###################################
export OMP_NUM_THREADS=1
#export EXE=$SLURM_SUBMIT_DIR/ns3d
#export EXE=$HOME/NS3D_active_brownian_bidisperse/trunk/ns3d
export EXE=$dirname/ns3d
#cd $SLURM_SUBMIT_DIR
echo '------------------------------------------------'
echo ' NODE USED = '$SLURM_NODELIST
echo ' SLURM_JOBID = '$SLURM_JOBID
echo ' OMP_NUM_THREADS = '$OMP_NUM_THREADS
echo ' CORES = '$cores
echo '------------------------------------------------'
ulimit -s unlimited
time srun -n $cores $EXE $WALLTIME >  ns3d_$cores\_$SLURM_JOBID.log
cd $dirname
ls -al
echo '------------------------------------------------'
echo ' COMPRESSING DATA IN : '$outname/$SLURM_JOBID.tar
echo '------------------------------------------------'
tar -cvzpf $SLURM_JOBID.tar *
mv $SLURM_JOBID.tar $outname


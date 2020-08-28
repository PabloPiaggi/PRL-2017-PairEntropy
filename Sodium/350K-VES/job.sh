#!/bin/bash
#SBATCH --ntasks=24             # total number of tasks across all nodes
#SBATCH --cpus-per-task=1       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=500M      # memory per cpu-core (4G is default)
#SBATCH --time=24:00:00         # total run time limit (HH:MM:SS)
#SBATCH --job-name="na-ves-350" 
#SBATCH --constraint=haswell|broadwell|skylake|cascade   # exclude ivy nodes

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PLUMED_NUM_THREADS=$SLURM_CPUS_PER_TASK

pwd; hostname; date

module purge
module load intel-mpi intel

############################################################################
# Variables definition
############################################################################
LAMMPS_HOME=/home/ppiaggi/Programs/Lammps/lammps-git-cpu/build3
LAMMPS_EXE=${LAMMPS_HOME}/lmp_della
cycles=1
threads_per_partition=4
############################################################################

############################################################################
# Run
############################################################################
if [ -e runno ] ; then
   #########################################################################
   # Restart runs
   #########################################################################
   nn=`tail -n 1 runno | awk '{print $1}'`
   srun $LAMMPS_EXE -partition 6x${threads_per_partition} -sf omp -sf intel -in Restart.lmp
   #########################################################################
else
   #########################################################################
   # First run
   #########################################################################
   nn=1
   # Number of partitions
   srun $LAMMPS_EXE -partition 6x${threads_per_partition} -sf omp -sf intel -in start.lmp
   #########################################################################
fi
############################################################################


############################################################################
# Prepare next run
############################################################################
# Back up
for j in $(seq 0 5)
do
        cp restart2.${j} restart2.${j}.${nn}
        cp restart.${j} restart.${j}.${nn}
        cp data.final.${j} data.final.${j}.${nn}
done
############################################################################

############################################################################
# Check number of cycles
############################################################################
mm=$((nn+1))
echo ${mm} > runno
#cheking number of cycles
if [ ${nn} -ge ${cycles} ]; then
  exit
fi
############################################################################

############################################################################
# Resubmitting again
############################################################################
sbatch < job.sh
############################################################################

date

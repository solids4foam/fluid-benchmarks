#!/bin/bash
#SBATCH --job-name=streamlines
#SBATCH --output=log.%x_%j.out
#SBATCH --error=log.%x_%j.err
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

# Good practice to navigate to the submission directory
cd $SLURM_SUBMIT_DIR

# Set number of threads for programs that use OpenMP parallelisation
export OMP_NUM_THREADS=$SLURM_NTASKS
#export OPENBLAS_NUM_THREADS=$SLURM_NTASKS
#export NUMEXPR_NUM_THREADS=$SLURM_NTASKS
#export MKL_NUM_THREADS=$SLURM_NTASKS
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Execute the run script
pvbatch --opengl-window-backend=EGL streamlines.py &> log.streamlines

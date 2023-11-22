#!/bin/bash
#

# Settings for batch job on UCL Legion:
#$ -l h_rt=2:00:00
#$ -l h_vmem=2G
#$ -N uptrop_test_job_source
#$ -pe smp 1
#$ -V -cwd
# -P hpc.10
# -l paid=1

# Load environment modules
#export OMP_NUM_THREADS=1

# Move to relevant working directory that has the Python script:
cd /lustre/scratch/scratch/ucfauxx/projects/uptrop/

# Activate your virtual environment:
source /home/ucfauxx/miniconda3/etc/profile.d/conda.sh 
conda activate uptrop
conda env list




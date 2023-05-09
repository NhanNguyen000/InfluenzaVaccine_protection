#!/bin/bash
#SBATCH --job-name=test1_2        # Name of job
#SBATCH --output=test1_2.out        # stdout
#SBATCH --error=test1_2.err         # stderr
#SBATCH --partition=cpu           # partition to use (check with sinfo)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks | Alternative: --ntasks-per-node
#SBATCH --threads-per-core=1      # Ensure we only get one logical CPU per core
#SBATCH --cpus-per-task=1         # Number of cores per task
#SBATCH --mem=16G                 # Memory per node | Alternative: --mem-per-cpu
#SBATCH --clusters=bioinf

cd /vol/projects/CIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen
Rscript 03_a07_responderScore_elasticModel_v2.R

/usr/bin/hostname


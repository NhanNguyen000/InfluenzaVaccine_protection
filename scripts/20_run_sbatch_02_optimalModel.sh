#!/bin/bash
#SBATCH --job-name=elasticNet_impVar       # Name of job
#SBATCH --output=scripts/elasticNet_optimal.out        # stdout
#SBATCH --error=scripts/elasticNet_optimal.err         # stderr
#SBATCH --partition=cpu           # partition to use (check with sinfo)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks | Alternative: --ntasks-per-node
#SBATCH --threads-per-core=1      # Ensure we only get one logical CPU per core
#SBATCH --cpus-per-task=1         # Number of cores per task
#SBATCH --mem=16G                 # Memory per node | Alternative: --mem-per-cpu
#SBATCH --time=30:00:00            # wall time limit (HH:MM:SS)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhan.nguyen@helmholtz-hzi.de
#SBATCH --clusters=bioinf

cd /vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/
Rscript scripts/15_elasticModel_03_paraRange_04_alpha.R

/usr/bin/hostname


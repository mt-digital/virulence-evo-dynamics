#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --output=VirEvoDyn-main.out
#SBATCH --partition=serc,normal

# Usage: scripts/slurm/main.sh MINORITY_HOMOPHILY MAJORITY_HOMOPHILY MIN_GROUP_FRAC
# 
# This script will write data to the data/main directory in the home project
# directory, which should be a symlink to a directory in SCRATCH.

module --force purge
module load devel
module load julia/1.9

julia scripts/run_trials.jl main --min_homophily=$1 --maj_homophily=$2 --min_group_frac=$3

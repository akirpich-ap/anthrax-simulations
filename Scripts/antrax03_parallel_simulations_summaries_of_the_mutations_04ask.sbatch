#!/bin/sh
#SBATCH --job-name=AZ_ANTR
#SBATCH --output=/ufrc/cummings/akirpich/antrax_simulation/Scripts/SLURM_LOGS/antrax03_parallel_simulations_summaries_of_the_mutations_04ask.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=akirpich@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=96:00:00
#SBATCH --account=epi 
#SBATCH --qos=epi

pwd;hostname;dateless 

module load R

# Running R using Rscript command
Rscript /ufrc/cummings/akirpich/antrax_simulation/R_Code/antrax03_parallel_simulations_summaries_of_the_mutations_04ask.R




#!/bin/sh
#SBATCH --job-name=04ant02
#SBATCH --output=/ufrc/cummings/akirpich/antrax_simulation/Scripts/SLURM_LOGS/antrax02_parallel_simulations_generic_all_settings_based_on_table_04ask_log_array%A_job%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=akirpich@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=96:00:00
#SBATCH --account=epi 
#SBATCH --qos=epi-b
#SBATCH --array=1-1000

pwd;hostname;dateless 

## Scanario ID: corresponding row from the parameters table.
SCENARIOID=2
## Setting up the absolute working directory path that is passed to the R script.
WORKDIR=/ufrc/cummings/akirpich/antrax_simulation
## Relative path for the parameters table in relation to WORKDIR that is passed to the script.
PARAMETER=Design/parameters_table_extended_updated_final_fixed.csv

module load R

# Running R using Rscript command
Rscript /ufrc/cummings/akirpich/antrax_simulation/R_Code/antrax02_parallel_simulations_generic_all_settings_based_on_table_04ask.R $SCENARIOID $WORKDIR $PARAMETER ${SLURM_ARRAY_TASK_ID}




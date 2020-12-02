// © Alexander Kirpich akirpich@gsu.edu 
// https://publichealth.gsu.edu/profile/alexander-kirpich/

Anthrax simulations for the manuscript “Alteration of exosporium surface oligosaccharides: Evidence of convergent patho-evolution in Bacillus anthracis”

The repository contains the following folders:

Data – the original anthrax sequences used in the simulation
Design – the folder contains the design file for the simulation scenarios
Documentation – the documentation folder
Plots – the intermediate plots produced within the code
R_Code – the R code used for simulation
R_Data – the folder in which the summaries from the individual runs from High Performance Computer (i.e. cluster or HPC) are combined.
R_Output – the folder contains the summary outputs that were used to produce parts of Figure 1 Including the quantiles for multiple simulations. 

R_Code folder files have the following goals:
antrax01_generate_parameters_for_simulations_05ask.R - The goal of this code is to define mutation rates a mutations in cholera bacteria over time.

antrax02_parallel_simulations_generic_all_settings_based_on_table_04ask.R - The goal of this code is to write a mutations in Anthrax bacteria over time.

antrax03_parallel_simulations_summaries_of_the_mutations_04ask.R - Here the parallel scenario i is considered and scenarioID is determined from the design table.

antrax04_parallel_simulations_summaries_of_the_mutations_summary_plots_02ask.R - This file is used to produce the simulation summaries across multiple scenarios that are used within the simulations.


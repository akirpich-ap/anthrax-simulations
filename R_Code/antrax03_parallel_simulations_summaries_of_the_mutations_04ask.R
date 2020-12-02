# // © Alexander Kirpich akirpich@gsu.edu
# This file is used to produce the simulation summaries across multiple scenarios that are used wihtin the simulations.
rm(list=ls(all=TRUE))

# library(Matrix) for blog-diagonal matrixes creation and other matrix manipulations.
library(Matrix)
# library(matrixStats) for rowMedians function.
library(matrixStats)
# library(MASS) for MVN simulations.
library(MASS)
# This package is requred to run in RScript mode rathen than interactive mode.
library(methods)


# Libraries to work with fasta files.
# For readDNAStringSet()
library(seqinr)
# For write.fasta()
library(Biostrings)



# Reading parameters from the command line.
# This piece of code is used to aquire parameters from command line/in batch mode via sbatch submission of HiperGator.
# args <- commandArgs(TRUE)
# print(args)
# current_scenario <- as.integer(args[1])
# work_directory_path <- args[2]
# parameter_table_relative_path <- args[3]
# fasta_file_with_data_relative_path <- args[4]
# current_simulation_dataset_id <- as.integer(args[5])


# Debugging step to run on local machine instead instead of the code right above used for HiPer Gator.
# current_scenario <- 1
# work_directory_path  <- "C:/Users/akirpich/Dropbox/GSU Research/2019 - Cummings Norris Antrax"
work_directory_path <- "/ufrc/cummings/akirpich/antrax_simulation/"
parameter_table_relative_path <- "Design/parameters_table_extended_updated_final_fixed.csv"


# Setting up the working directory.
# This path is passed via work_directory_path variable that is passed via command line parameter on HiperGator.
setwd(work_directory_path)
# Extra check
getwd()



# Part 01
# Reading the table of settings inside R.
table_of_settings  <- read.table( parameter_table_relative_path,   sep = ",", header = TRUE )

# Fix 2019.04.16.
# Looping over scenario_id-s and over the realizations within each scenario.

# Generating additional variables for that.
# Total number of scenarios.
no_of_scenarios <- dim(table_of_settings)[1]
# Total number of replications within each scenario.
scenario_id <- 1
no_of_replications <- table_of_settings$no_of_replications[scenario_id]
# Number of samples in each replication and each scenario i.e. population size.
no_of_samples <- table_of_settings$parameter_n[scenario_id]


# Debuggins step 
# Defining smaller frame then the that is needed in case RAM is an issue. 2.2GB is needed for the complete frame.
# no_of_scenarios <- 1



# Fix 2019.09.30 
# Generating summaries table for indicators for presence of stop codons.
# beginning
stopcodon_indicators_summary_array_beginning <- array(0, dim = c(no_of_scenarios, no_of_replications, no_of_samples))
dim(stopcodon_indicators_summary_array_beginning)
# end
stopcodon_indicators_summary_array_end       <- array(0, dim = c(no_of_scenarios, no_of_replications, no_of_samples))
dim(stopcodon_indicators_summary_array_end)
# ever
stopcodon_indicators_summary_array_ever      <- array(0, dim = c(no_of_scenarios, no_of_replications, no_of_samples))
dim(stopcodon_indicators_summary_array_ever)



# Fix 2019.09.30 
# Generating summaries table for hamming distances.
# beginning
hamming_distance_summary_array_beginning <- array(0, dim = c(no_of_scenarios, no_of_replications, no_of_samples))
dim(hamming_distance_summary_array_beginning)
# end
hamming_distance_summary_array_end       <- array(0, dim = c(no_of_scenarios, no_of_replications, no_of_samples))
dim(hamming_distance_summary_array_end)



# Fix 2019.09.30 
# Generating summaries table for mutant/wildtype ratio.
# beginning
mutant_wildtype_summary_array_beginning <- array(0, dim = c(no_of_scenarios, no_of_replications, no_of_samples))
dim(mutant_wildtype_summary_array_beginning)
# end
mutant_wildtype_summary_array_end       <- array(0, dim = c(no_of_scenarios, no_of_replications, no_of_samples))
dim(mutant_wildtype_summary_array_end)





# Looping over scenarios
for ( scenario_id in c(1:no_of_scenarios) )
{
  # Debugging step
  # scenario_id <- 1

  
  for ( current_replication in c(1:table_of_settings$no_of_replications[scenario_id]  )  )
  {
    # Debugging step
    # current_replication <- 1

    
    # Stop codons
    
    # Generating current path for the stop codon presence over time status table
    stop_codon_presence_over_time_short_current_path   <- paste("R_Output/stop_codon_presence_over_time_short_scenario_id_", table_of_settings$scenario_id[scenario_id], "_replication_", current_replication, ".csv", sep ="")
    # Reading the stop codon presence over time status table
    stop_codon_presence_over_time_short_current   <- read.table( stop_codon_presence_over_time_short_current_path,   sep = ",", header = TRUE )
    # Debugging step
    # dim(stop_codon_presence_over_time_short_current)
    
    # Debugging step
    cat("  scenario_id = ", table_of_settings$scenario_id[scenario_id], " and replication = ", current_replication, "Stop Codon Status Was Read ! \n")
    
    # Checking for stop codon presence anywhere in the simulation
    stop_codon_presence_over_time_short_current_indicators  <-  (stop_codon_presence_over_time_short_current == "P")
    # Debubbing step
    # dim(stop_codon_presence_over_time_short_current_indicators) 

    # Extracting indicators
    # Extracting value for the first time within the simulation.
    stop_codon_presence_over_time_short_current_indicator_first <-  ( stop_codon_presence_over_time_short_current_indicators[, 1 ]  > 0) * 1
    # Extracting value for the  last time within the simulation.
    stop_codon_presence_over_time_short_current_indicator_last  <-  ( stop_codon_presence_over_time_short_current_indicators[, dim(stop_codon_presence_over_time_short_current_indicators)[2] ]  > 0) * 1
    # Extracting value when it ever appeared within the simulation.
    stop_codon_presence_over_time_short_current_indicator_ever  <-  ( rowSums(stop_codon_presence_over_time_short_current_indicators)  > 0) * 1
    
    # Debubbing step
    # length(indicators_for_current_simulation)
    
    # Saving the results into the data frame that was created earlier.
    # Fist time slots array
    dim(stopcodon_indicators_summary_array_beginning[scenario_id, , ])
    stopcodon_indicators_summary_array_beginning[scenario_id, current_replication, ] <- stop_codon_presence_over_time_short_current_indicator_first
    # Last time slots array
    dim(stopcodon_indicators_summary_array_end[scenario_id, , ])
    stopcodon_indicators_summary_array_end[scenario_id, current_replication, ]       <- stop_codon_presence_over_time_short_current_indicator_last
    # Ever time slots array
    dim(stopcodon_indicators_summary_array_ever[scenario_id, , ])
    stopcodon_indicators_summary_array_ever[scenario_id, current_replication, ]      <- stop_codon_presence_over_time_short_current_indicator_ever
    

    
    
    # Hamming Distance
    
    # Generating current path for the stop codon presence over time status table
    hamming_distance_over_time_short_current_path   <- paste("R_Output/hamming_distance_over_time_short_scenario_id_", table_of_settings$scenario_id[scenario_id], "_replication_", current_replication, "_short.csv", sep ="")
    # Reading the stop codon presence over time status table
    hamming_distance_over_time_short_current   <- read.table( hamming_distance_over_time_short_current_path,   sep = ",", header = TRUE )
    # Debugging step
    # dim(hamming_distance_over_time_short_current)
    
    # Debugging step
    cat("  scenario_id = ", table_of_settings$scenario_id[scenario_id], " and replication = ", current_replication, "Hamming Distance Was Read ! \n")
    
    # Extracting values
    # Extracting value for the first time within the simulation.
    hamming_distance_over_time_short_current_indicator_first <-  hamming_distance_over_time_short_current[, 1 ] 
    # Extracting value for the  last time within the simulation.
    hamming_distance_over_time_short_current_indicator_last  <-  hamming_distance_over_time_short_current[, dim(hamming_distance_over_time_short_current)[2] ]

    # Debubbing step
    # length(hamming_distance_over_time_short_current)
    
    # Saving the results into the data frame that was created earlier.
    # Fist time slots array
    dim(hamming_distance_summary_array_beginning[scenario_id, , ])
    hamming_distance_summary_array_beginning[scenario_id, current_replication, ] <- hamming_distance_over_time_short_current_indicator_first
    # Last time slots array
    dim(hamming_distance_summary_array_end[scenario_id, , ])
    hamming_distance_summary_array_end[scenario_id, current_replication, ]       <- hamming_distance_over_time_short_current_indicator_last

    

    
    # Mutant/Wildtype
    
    # Generating current path for the stop codon presence over time status table
    wildtype_mutant_status_over_time_short_current_path   <- paste("R_Output/wildtype_mutant_status_over_time_short_scenario_id_", table_of_settings$scenario_id[scenario_id], "_replication_", current_replication, "_short.csv", sep ="")
    # Reading the stop codon presence over time status table
    wildtype_mutant_status_over_time_short_current   <- read.table( wildtype_mutant_status_over_time_short_current_path,   sep = ",", header = TRUE )
    # Debugging step
    # dim(wildtype_mutant_status_over_time_short_current)
    
    # Debugging step
    cat("  scenario_id = ", table_of_settings$scenario_id[scenario_id], " and replication = ", current_replication, "Mutant/Wildtype Status Was Read ! \n")
    
    # Checking for stop codon presence anywhere in the simulation
    wildtype_mutant_status_over_time_short_current_indicators  <-  (wildtype_mutant_status_over_time_short_current == "M")
    # Debubbing step
    # dim(wildtype_mutant_status_over_time_short_current_indicators) 
    
    # Extracting indicators
    # Extracting value for the first time within the simulation.
    wildtype_mutant_status_over_time_short_current_indicator_first <-  ( wildtype_mutant_status_over_time_short_current_indicators[, 1 ]  > 0) * 1
    # Extracting value for the  last time within the simulation.
    wildtype_mutant_status_over_time_short_current_indicator_last  <-  ( wildtype_mutant_status_over_time_short_current_indicators[, dim(wildtype_mutant_status_over_time_short_current_indicators)[2] ]  > 0) * 1


    # Saving the results into the data frame that was created earlier.
    # Fist time slots array
    dim(mutant_wildtype_summary_array_beginning[scenario_id, , ])
    mutant_wildtype_summary_array_beginning[scenario_id, current_replication, ] <- wildtype_mutant_status_over_time_short_current_indicator_first
    # Last time slots array
    dim(mutant_wildtype_summary_array_end[scenario_id, , ])
    mutant_wildtype_summary_array_end[scenario_id, current_replication, ]       <- wildtype_mutant_status_over_time_short_current_indicator_last


        
  # End of -> for (scenario_id in c(1:no_of_scenarios) )
  }     
  
# End of -> for ( scenario_id in c(1:no_of_scenarios) )  
}  





# Fix 2019.11.07.
# Exporting the large results.

# Stop codons
# beginning
stopcodon_indicators_summary_array_beginning_path <- paste("R_Data/stopcodon_indicators_summary_array_beginning.RData")
save( stopcodon_indicators_summary_array_beginning,  file = stopcodon_indicators_summary_array_beginning_path )
# end
stopcodon_indicators_summary_array_end_path       <- paste("R_Data/stopcodon_indicators_summary_array_end.RData")
save( stopcodon_indicators_summary_array_end,        file = stopcodon_indicators_summary_array_end_path )
# ever
stopcodon_indicators_summary_array_ever_path      <- paste("R_Data/stopcodon_indicators_summary_array_ever.RData")
save( stopcodon_indicators_summary_array_ever,       file = stopcodon_indicators_summary_array_ever_path )


# Hamming distance
# beginning
hamming_distance_summary_array_beginning_path <- paste("R_Data/hamming_distance_summary_array_beginning.RData")
save( hamming_distance_summary_array_beginning,  file = hamming_distance_summary_array_beginning_path )
# end
hamming_distance_summary_array_end_path       <- paste("R_Data/hamming_distance_summary_array_end.RData")
save( hamming_distance_summary_array_end,        file = hamming_distance_summary_array_end_path )


# Mutant/Wildtype
# beginning
mutant_wildtype_summary_array_beginning_path <- paste("R_Data/mutant_wildtype_summary_array_beginning.RData")
save( mutant_wildtype_summary_array_beginning,  file = mutant_wildtype_summary_array_beginning_path )
# end
mutant_wildtype_summary_array_end_path       <- paste("R_Data/mutant_wildtype_summary_array_end.RData")
save( mutant_wildtype_summary_array_end,        file = mutant_wildtype_summary_array_end_path )










# Generating vector for probabilities.
probability_quantiles <- c(0, 0.024, 0.25, 0.50, 0.75, 0.976, 1)

# Generating summaries table for quantiles across all runs.

# Stop codons
# beginning
stop_codon_summary_table_beginning <- data.frame( matrix(0, nrow = no_of_scenarios, ncol = length(probability_quantiles) + 1 ) )
names(stop_codon_summary_table_beginning) <- c("Type", "0%", "2.4%", "25%", "50%", "75%", "97.6%", "100%" )
dim(stop_codon_summary_table_beginning)
# end
stop_codon_summary_table_end <- data.frame( matrix(0, nrow = no_of_scenarios, ncol = length(probability_quantiles) + 1 ) )
names(stop_codon_summary_table_end) <- c("Type", "0%", "2.4%", "25%", "50%", "75%", "97.6%", "100%" )
dim(stop_codon_summary_table_end)
# ever
stop_codon_summary_table_ever <- data.frame( matrix(0, nrow = no_of_scenarios, ncol = length(probability_quantiles) + 1 ) )
names(stop_codon_summary_table_ever) <- c("Type", "0%", "2.4%", "25%", "50%", "75%", "97.6%", "100%" )
dim(stop_codon_summary_table_ever)


# hamming distance
# beginning
hamming_distance_summary_table_beginning <- data.frame( matrix(0, nrow = no_of_scenarios, ncol = length(probability_quantiles) + 1 ) )
names(hamming_distance_summary_table_beginning) <- c("Type", "0%", "2.4%", "25%", "50%", "75%", "97.6%", "100%" )
dim(hamming_distance_summary_table_beginning)
# end
hamming_distance_summary_table_end <- data.frame( matrix(0, nrow = no_of_scenarios, ncol = length(probability_quantiles) + 1 ) )
names(hamming_distance_summary_table_end) <- c("Type", "0%", "2.4%", "25%", "50%", "75%", "97.6%", "100%" )
dim(hamming_distance_summary_table_end)


# mutant-wildtype
# beginning
mutant_wildtype_summary_table_beginning <- data.frame( matrix(0, nrow = no_of_scenarios, ncol = length(probability_quantiles) + 1 ) )
names(mutant_wildtype_summary_table_beginning) <- c("Type", "0%", "2.4%", "25%", "50%", "75%", "97.6%", "100%" )
dim(mutant_wildtype_summary_table_beginning)
# end
mutant_wildtype_summary_table_end <- data.frame( matrix(0, nrow = no_of_scenarios, ncol = length(probability_quantiles) + 1 ) )
names(mutant_wildtype_summary_table_end) <- c("Type", "0%", "2.4%", "25%", "50%", "75%", "97.6%", "100%" )
dim(mutant_wildtype_summary_table_end)







# Fix 2019.10.01
# Doing two summaries simulteneously: the first summary if for proportion quantiles and the second summary is for the median hamming distance quantiles.
# Looping over scenarios
for ( scenario_id in c(1:no_of_scenarios) )
{
  # Debugging step
  # scenario_id <- 1
  
  # Computing summaries and qunatiles.
  
   
  
  # Stop codons
  
  # beginning
  # Extracting summaries for the current scenario.
  stopcodon_indicators_summary_array_beginning_scenario_id <- rowSums( stopcodon_indicators_summary_array_beginning[scenario_id, , ] ) / dim(stopcodon_indicators_summary_array_beginning[scenario_id, , ])[2]
  # length(stopcodon_indicators_summary_array_beginning_scenario_id)
  # Saving the results into the table.
  stop_codon_summary_table_beginning[scenario_id, c(2:dim(stop_codon_summary_table_beginning)[2]) ] <- quantile(  stopcodon_indicators_summary_array_beginning_scenario_id,  probs = probability_quantiles  )  
  stop_codon_summary_table_beginning[scenario_id, 1 ] <-  as.character(table_of_settings$scenario_id[scenario_id])

  # end
  # Extracting summaries for the current scenario.
  stopcodon_indicators_summary_array_end_scenario_id <- rowSums( stopcodon_indicators_summary_array_end[scenario_id, , ] ) / dim(stopcodon_indicators_summary_array_end[scenario_id, , ])[2]
  # length(stopcodon_indicators_summary_array_end_scenario_id)
  # Saving the results into the table.
  stop_codon_summary_table_end[scenario_id, c(2:dim(stop_codon_summary_table_end)[2]) ] <- quantile(  stopcodon_indicators_summary_array_end_scenario_id,  probs = probability_quantiles  )  
  stop_codon_summary_table_end[scenario_id, 1 ] <-  as.character(table_of_settings$scenario_id[scenario_id])

  # ever
  # Extracting summaries for the current scenario.
  stopcodon_indicators_summary_array_ever_scenario_id <- rowSums( stopcodon_indicators_summary_array_ever[scenario_id, , ] ) / dim(stopcodon_indicators_summary_array_ever[scenario_id, , ])[2]
  # length(stopcodon_indicators_summary_array_ever_scenario_id)
  # Saving the results into the table.
  stop_codon_summary_table_ever[scenario_id, c(2:dim(stop_codon_summary_table_ever)[2]) ] <- quantile(  stopcodon_indicators_summary_array_ever_scenario_id,  probs = probability_quantiles  )  
  stop_codon_summary_table_ever[scenario_id, 1 ] <-  as.character(table_of_settings$scenario_id[scenario_id])
  
  
  
  # Hamming distance
  
  # beginning
  # Extracting summaries for the current scenario.
  hamming_distance_summary_array_beginning_scenario_id <- rowSums( hamming_distance_summary_array_beginning[scenario_id, , ] ) / dim(hamming_distance_summary_array_beginning[scenario_id, , ])[2]
  # length(hamming_distance_summary_array_beginning_scenario_id)
  # Saving the results into the table.
  hamming_distance_summary_table_beginning[scenario_id, c(2:dim(hamming_distance_summary_table_beginning)[2]) ] <- quantile(  hamming_distance_summary_array_beginning_scenario_id,  probs = probability_quantiles  )  
  hamming_distance_summary_table_beginning[scenario_id, 1 ] <-  as.character(table_of_settings$scenario_id[scenario_id])
  
  # end
  # Extracting summaries for the current scenario.
  hamming_distance_summary_array_end_scenario_id <- rowSums( hamming_distance_summary_array_end[scenario_id, , ] ) / dim(hamming_distance_summary_array_end[scenario_id, , ])[2]
  # length(hamming_distance_summary_array_end_scenario_id)
  # Saving the results into the table.
  hamming_distance_summary_table_end[scenario_id, c(2:dim(hamming_distance_summary_table_end)[2]) ] <- quantile(  hamming_distance_summary_array_end_scenario_id,  probs = probability_quantiles  )  
  hamming_distance_summary_table_end[scenario_id, 1 ] <-  as.character(table_of_settings$scenario_id[scenario_id])
  
  
  
  # Mutant/Wildtype
  
  # beginning
  # Extracting summaries for the current scenario.
  mutant_wildtype_summary_array_beginning_scenario_id <- rowSums( mutant_wildtype_summary_array_beginning[scenario_id, , ] ) / dim(mutant_wildtype_summary_array_beginning[scenario_id, , ])[2]
  # length(mutant_wildtype_summary_array_beginning_scenario_id)
  # Saving the results into the table.
  mutant_wildtype_summary_table_beginning[scenario_id, c(2:dim(mutant_wildtype_summary_table_beginning)[2]) ] <- quantile(  mutant_wildtype_summary_array_beginning_scenario_id,  probs = probability_quantiles  )  
  mutant_wildtype_summary_table_beginning[scenario_id, 1 ] <-  as.character(table_of_settings$scenario_id[scenario_id])
  
  # end
  # Extracting summaries for the current scenario.
  mutant_wildtype_summary_array_end_scenario_id <- rowSums( mutant_wildtype_summary_array_end[scenario_id, , ] ) / dim(mutant_wildtype_summary_array_end[scenario_id, , ])[2]
  # length(mutant_wildtype_summary_array_end_scenario_id)
  # Saving the results into the table.
  mutant_wildtype_summary_table_end[scenario_id, c(2:dim(mutant_wildtype_summary_table_end)[2]) ] <- quantile(  mutant_wildtype_summary_array_end_scenario_id,  probs = probability_quantiles  )  
  mutant_wildtype_summary_table_end[scenario_id, 1 ] <-  as.character(table_of_settings$scenario_id[scenario_id])
  

# End of -> for ( scenario_id in c(1:no_of_scenarios) )
}  




# Exporting the results.

# Stop codons
# beginning
stop_codon_summary_table_beginning_path <- paste("R_Output/stop_codon_summary_table_beginning.csv")
write.csv( x = stop_codon_summary_table_beginning,  file = stop_codon_summary_table_beginning_path, row.names = FALSE, quote = FALSE )
# end
stop_codon_summary_table_end_path       <- paste("R_Output/stop_codon_summary_table_end.csv")
write.csv( x = stop_codon_summary_table_end,        file = stop_codon_summary_table_end_path, row.names = FALSE, quote = FALSE )
# ever
stop_codon_summary_table_ever_path      <- paste("R_Output/stop_codon_summary_table_ever.csv")
write.csv( x = stop_codon_summary_table_ever,       file = stop_codon_summary_table_ever_path, row.names = FALSE, quote = FALSE )


# Hamming distance
# beginning
hamming_distance_summary_table_beginning_path <- paste("R_Output/hamming_distance_summary_table_beginning.csv")
write.csv( x = hamming_distance_summary_table_beginning,  file = hamming_distance_summary_table_beginning_path, row.names = FALSE, quote = FALSE )
# end
hamming_distance_summary_table_end_path       <- paste("R_Output/hamming_distance_summary_table_end.csv")
write.csv( x = hamming_distance_summary_table_end,        file = hamming_distance_summary_table_end_path, row.names = FALSE, quote = FALSE )


# Mutant/Wildtype
# beginning
mutant_wildtype_summary_table_beginning_path <- paste("R_Output/mutant_wildtype_summary_table_beginning.csv")
write.csv( x = mutant_wildtype_summary_table_beginning,  file = mutant_wildtype_summary_table_beginning_path, row.names = FALSE, quote = FALSE )
# end
mutant_wildtype_summary_table_end_path       <- paste("R_Output/mutant_wildtype_summary_table_end.csv")
write.csv( x = mutant_wildtype_summary_table_end,        file = mutant_wildtype_summary_table_end_path, row.names = FALSE, quote = FALSE )










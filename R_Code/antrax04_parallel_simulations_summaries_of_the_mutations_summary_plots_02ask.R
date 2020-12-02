# // © Alexander Kirpich akirpich@gsu.edu
# This code is used for addiitonal plotting.
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
# For integer.frequency(x,bins)
library(prettyR)


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
work_directory_path  <- "C:/Users/akirpich/Dropbox/GSU Research/2019 - Cummings Norris Antrax"
parameter_table_relative_path <- "Design/parameters_table_extended_updated_final_fixed.csv"


# Setting up the working directory.
# This path is passed via work_directory_path variable that is passed via command line parameter on HiperGator.
setwd(work_directory_path)
# Extra check
getwd()

# Reading the table of settings inside R.
table_of_settings  <- read.table( parameter_table_relative_path,   sep = ",", header = TRUE )

# Generating additional variables for that.
# Total number of scenarios.
no_of_scenarios <- dim(table_of_settings)[1]
# Total number of replications within each scenario.
# scenario_id <- 1
no_of_replications <- table_of_settings$no_of_replications[scenario_id]
# Number of samples in each replication and each scenario i.e. population size.
no_of_samples <- table_of_settings$parameter_n[scenario_id]





# Importing the results

# Hamming distance
# beginning
hamming_distance_summary_table_beginning_path <- paste("R_Output/hamming_distance_summary_table_beginning.csv")
hamming_distance_summary_table_beginning      <- read.table(  file = hamming_distance_summary_table_beginning_path, sep = "," , header = TRUE )
# end
hamming_distance_summary_table_end_path       <- paste("R_Output/hamming_distance_summary_table_end.csv")
hamming_distance_summary_table_end            <- read.table(  file = hamming_distance_summary_table_end_path, sep = "," , header = TRUE )

# Fix 2019.11.12
# Hamming distance distribution
# beginning
hamming_distance_summary_array_beginning_path <- paste("R_Data/hamming_distance_summary_array_beginning.RData")
load( hamming_distance_summary_array_beginning_path )
# end
hamming_distance_summary_array_end_path       <- paste("R_Data/hamming_distance_summary_array_end.RData")
load( hamming_distance_summary_array_end_path )





# Mutant/Wildtype
# beginning
mutant_wildtype_summary_table_beginning_path <- paste("R_Output/mutant_wildtype_summary_table_beginning.csv")
mutant_wildtype_summary_table_beginning      <- read.table(  file = mutant_wildtype_summary_table_beginning_path, sep = "," , header = TRUE )
# end
mutant_wildtype_summary_table_end_path <- paste("R_Output/mutant_wildtype_summary_table_end.csv")
mutant_wildtype_summary_table_end      <- read.table(  file = mutant_wildtype_summary_table_end_path, sep = "," , header = TRUE )




# Mutant-Wildtype Proportions

# Extracting data
mutant_wildtype_summary_table_end_median   <- mutant_wildtype_summary_table_end[, c(1,6) ]
mutant_wildtype_summary_table_end_low_high <- mutant_wildtype_summary_table_end[, c(1,3,6,8) ]


t_mutant_wildtype_summary_table_end_median      <- t(mutant_wildtype_summary_table_end_median)
t_mutant_wildtype_summary_table_end_median_low  <- rbind(    as.numeric( t(mutant_wildtype_summary_table_end_median)[2,] ), 
                                                         1 - as.numeric( t(mutant_wildtype_summary_table_end_median)[2,] ) )
colnames(t_mutant_wildtype_summary_table_end_median_low) <- t_mutant_wildtype_summary_table_end_median[1,]
# rownames(t_mutant_wildtype_summary_table_end_median_low) <- "Proportion"


# Producing plot for the "current_parameter" splittede b/w treatment and control groups.
pdf( paste("Plots/antrax04_parallel_simulations_summaries_of_the_mutations_summary_plots_01ask_proportions.pdf", sep=""), height = 6 , width = 8 )
barplot( t_mutant_wildtype_summary_table_end_median_low, col=c("orange", "darkblue"), 
         legend = TRUE, border =  TRUE, xlim = c(0,15), args.legend = list(bty="n", border=TRUE), 
         ylab = "Proportion of Bacteria With Triplets", xlab="Scenario ID",
         names.arg = colnames(t_mutant_wildtype_summary_table_end_median_low), cex.names=0.75 )
legend( "topright", inset= c(0, 0.4), legend = c("Mutated", "Not Mutated"), col = "black", fill = c("orange", "darkblue"),   pt.cex = 2,  cex = 0.95  )  
dev.off()





# Hamming distance

# Extracting data
hamming_distance_summary_table_end_median   <- hamming_distance_summary_table_end[, c(1,6) ]
hamming_distance_summary_table_end_low_high <- hamming_distance_summary_table_end[, c(1,3,6,8) ]


t_hamming_distance_summary_table_end_median      <- t(hamming_distance_summary_table_end_median)
t_hamming_distance_summary_table_end_median_low  <- rbind(          as.numeric( t(hamming_distance_summary_table_end_median)[2,] ), 
                                                             2000 - as.numeric( t(hamming_distance_summary_table_end_median)[2,] ) )
colnames(t_hamming_distance_summary_table_end_median_low) <- t_hamming_distance_summary_table_end_median[1,]
# rownames(t_hamming_distance_summary_table_end_median_low) <- "Proportion"


# Producing plot for the "current_parameter" splittede b/w treatment and control groups.
pdf( paste("Plots/antrax04_parallel_simulations_summaries_of_the_mutations_summary_plots_01ask_counts.pdf", sep=""), height = 6 , width = 8 )
barplot( t_hamming_distance_summary_table_end_median_low, col=c("orange", "darkblue"), 
         legend = TRUE, border =  TRUE, xlim = c(0,15), args.legend = list(bty="n", border=TRUE), 
         ylab = "Average Counts that are Different", xlab="Scenario ID",
         names.arg = colnames(t_hamming_distance_summary_table_end_median_low), cex.names=0.75 )
legend( "topright", inset= c(0, 0.4), legend = c("Mutated", "Not Mutated"), col = "black", fill = c("orange", "darkblue"),   pt.cex = 2,  cex = 0.95  )  
dev.off()



no_of_time_slots <- 100
parameter_interval_time_mean <- 7

list_of_active_years <- seq( from = 1, to = no_of_time_slots, by = parameter_interval_time_mean ) 

x_list = c( 1:no_of_time_slots )
y_list = rep(0.1, no_of_time_slots )
y_list[ list_of_active_years ] <- 1

# Producing plot for the "current_parameter" splittede b/w treatment and control groups.
pdf( paste("Plots/antrax04_parallel_simulations_summaries_of_the_mutations_summary_plots_01ask_spikes.pdf", sep=""), height = 6 , width = 16 )
plot( x = x_list, y = y_list, col='blue4', lwd=3, type = "l",
      xlab= "Time",
      ylab= "Probability Spikes Example",
      yaxt="n")
dev.off()











# Fix 2019.11.12.
# Summaries for Dr. Norris upon request.
# Getting the maximum and minimum possible Hamming distances across all runs.
max_hamming_distance <- max( hamming_distance_summary_array_end )
min_hamming_distance <- min( hamming_distance_summary_array_end )
# Generating a range
range_hamming_distance <- c( min_hamming_distance: max_hamming_distance )


# Creating summary table
hamming_distance_summary_end_matrix <- matrix(0, nrow = no_of_scenarios, ncol = length(range_hamming_distance) )
dim(hamming_distance_summary_end_matrix)
colnames(hamming_distance_summary_end_matrix) <- range_hamming_distance
rownames(hamming_distance_summary_end_matrix) <- table_of_settings$scenario_id






# Looping over scenarios
for ( scenario_id in c(1:no_of_scenarios) )
{
  # Debugging step
  # scenario_id <- 1
  # scenario_id <- 10
  
  # Hamming Distance
  
  # Etracting the current scenario table
  hamming_distance_summary_current_scenario <- hamming_distance_summary_array_end[ scenario_id, , ]
  
  # Summarizing the average over 10000 replications
  current_averaged_summary_table <- table(hamming_distance_summary_current_scenario)/table_of_settings$no_of_replications[scenario_id]
  
  # Debuggins step
  cat("Scenario ID ", scenario_id, "\n")
  
  # matching names 
  available_names <- names( current_averaged_summary_table )
  available_names_integers <- as.integer( available_names )
  
  # Debuggins step
  cat("available_names_integers ", available_names_integers, "\n")
  
  # Saving the results
  hamming_distance_summary_end_matrix[scenario_id, available_names_integers+1] <- current_averaged_summary_table
  
  # Debuggins step
  cat("current_averaged_summary_table ", current_averaged_summary_table, "\n")
  
  # Debuggins step
  cat(   "hamming_distance_summary_end_matrix[scenario_id, ]\n" )
  print(  hamming_distance_summary_end_matrix[scenario_id, ] )

    

  # End of -> for ( scenario_id in c(1:no_of_scenarios) )  
}  





hamming_distance_summary_end_matrix_percent <- hamming_distance_summary_end_matrix/100

# Exporing the results 
# Regular
hamming_distance_summary_end_matrix_percent_text <- cbind( rownames(hamming_distance_summary_end_matrix_percent),  hamming_distance_summary_end_matrix_percent )
colnames(hamming_distance_summary_end_matrix_percent_text)[1]  <- "Type"
hamming_distance_summary_end_matrix_percent_text_path <- paste("R_Output/hamming_distance_summary_end_matrix_percent_text.csv")
write.csv( x = hamming_distance_summary_end_matrix_percent_text,  file = hamming_distance_summary_end_matrix_percent_text_path, row.names = FALSE, quote = FALSE )
# Percent
hamming_distance_summary_end_matrix_text <- cbind( rownames(hamming_distance_summary_end_matrix),  hamming_distance_summary_end_matrix )
colnames(hamming_distance_summary_end_matrix_text)[1]  <- "Type"
hamming_distance_summary_end_matrix_text_path <- paste("R_Output/hamming_distance_summary_end_matrix_text.csv")
write.csv( x = hamming_distance_summary_end_matrix_text,  file = hamming_distance_summary_end_matrix_text_path, row.names = FALSE, quote = FALSE )





# Generating plots for different scenarios
pdf( paste("Plots/antrax04_parallel_simulations_summaries_of_the_mutations_summary_plots_01ask_percent_points.pdf", sep=""), height = 16 , width = 10 )
par( mfrow = c( no_of_scenarios/2, 2),  mar=c(6.2, 4.5, 2, 2), oma=c(2, 1.7, 5, 2)  )
# par( mfrow = c( no_of_scenarios/2, 2)  )

# Looping over scenarios
for ( scenario_id in c(1:no_of_scenarios) )
{
  # Plotting the corresponding parameter values
  plot( x = colnames(hamming_distance_summary_end_matrix_percent),
        y = hamming_distance_summary_end_matrix_percent[scenario_id,],
        # Plot limits.
        xlim = c(min(as.integer(colnames(hamming_distance_summary_end_matrix_percent))), max(as.integer(colnames(hamming_distance_summary_end_matrix_percent))) ),
        #ylim = c(0, max(hamming_distance_summary_end_matrix_percent)* 1.03),
        ylim = c(-1 , max(hamming_distance_summary_end_matrix_percent[scenario_id,])* 1.03),
        # separation by color.
        # type = "l",
        col  = "blue4",
        # Points width and type.
        lwd = 3,
        pch = 19,
        cex = 2,
        # Main title and axis titles
        main = paste(table_of_settings$scenario_id[scenario_id]),
        xlab = "Distance",
        ylab = "Percent",
        # removing labels from y=axies
        #xaxt='n',
        cex.axis = 1.7,
        cex.lab  = 1.7,
        cex.main = 1.7
  )
 
# End of -> for ( scenario_id in c(1:no_of_scenarios) )   
}
title("Distribution of Hamming Distances Averages From 1000 Replications\n Across Bacteria Population (in Percent)", outer=TRUE, cex.main = 2) 
dev.off()



# Generating plots for different scenarios
pdf( paste("Plots/antrax04_parallel_simulations_summaries_of_the_mutations_summary_plots_01ask_percent_lines.pdf", sep=""), height = 16 , width = 10 )
par( mfrow = c( no_of_scenarios/2, 2),  mar=c(6.2, 4.5, 2, 2), oma=c(2, 1.7, 5, 2)  )
# par( mfrow = c( no_of_scenarios/2, 2)  )

# Looping over scenarios
for ( scenario_id in c(1:no_of_scenarios) )
{
  # Plotting the corresponding parameter values
  plot( x = colnames(hamming_distance_summary_end_matrix_percent),
        y = hamming_distance_summary_end_matrix_percent[scenario_id,],
        # Plot limits.
        xlim = c(min(as.integer(colnames(hamming_distance_summary_end_matrix_percent))), max(as.integer(colnames(hamming_distance_summary_end_matrix_percent))) ),
        #ylim = c(0, max(hamming_distance_summary_end_matrix_percent)* 1.03),
        ylim = c(-1 , max(hamming_distance_summary_end_matrix_percent[scenario_id,])* 1.03),
        # separation by color.
        type = "l",
        col  = "blue4",
        # Points width and type.
        lwd = 3,
        pch = 19,
        cex = 2,
        # Main title and axis titles
        main = paste(table_of_settings$scenario_id[scenario_id]),
        xlab = "Distance",
        ylab = "Percent",
        # removing labels from y=axies
        #xaxt='n',
        cex.axis = 1.7,
        cex.lab  = 1.7,
        cex.main = 1.7
  )
  
  # End of -> for ( scenario_id in c(1:no_of_scenarios) )   
}
title("Distribution of Hamming Distances Averages From 1000 Replications\n Across Bacteria Population (in Percent)", outer=TRUE, cex.main = 2) 
dev.off()

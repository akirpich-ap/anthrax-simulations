# // © Alexander Kirpich akirpich@gsu.edu
# The goal of this code is to define mutation rates a mutations in cholera bacteria over time.
rm(list=ls(all=TRUE))

# 20 Digits Precision Representation
options(scipen=20)


# library(Matrix) for blog-diagonal matrixes creation and other matrix manipulations.
library(Matrix)
# library(MASS) for MVN simulations.
library(MASS)
# This package is requred to run in RScript mode rathen than interactive mode.
library(methods)


# Libraries to work with fasta files.
# For readDNAStringSet()
library(seqinr)
# For write.fasta()
library(Biostrings)


# Loading package requred to read library(readxl)
library(readxl)


# Debugging step to run on local machine instead instead of the code right above used for HiPer Gator.
# work_directory_path  <- "E:/Alexander/UF Research/2019 - Cummings Norris Antrax"
work_directory_path  <- "C:/Users/akirpich/Dropbox/GSU Research/2019 - Cummings Norris Antrax"

# Setting up the working directory.
setwd(work_directory_path)
# Extra check
getwd()


parameter_table_relative_path <- "Design/parameters_table_extended_updated_final_fixed.csv"
parameter_table <- read.csv(parameter_table_relative_path)
print(parameter_table)
names(parameter_table)


# Adding source fasta file path
parameter_table$source_fasta_path <- rep( "Data/anthoseoperonWTshort", dim(parameter_table)[1] )



# Fix 2019.09.11.
# Fixing the multipliers for the gentimeslot
# We are using the paper to define the number of generations per year.
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000461
# sporadic
# Assume just 1 generation for now
parameter_table$parameter_gentimeslot_wildtype_sporadic <- rep(1, dim(parameter_table)[1] )
parameter_table$parameter_gentimeslot_mutant_sporadic   <- rep(1, dim(parameter_table)[1] )
# Active
# Assume 43 based on the paper.
parameter_table$parameter_gentimeslot_wildtype_active   <- rep(43, dim(parameter_table)[1] )
parameter_table$parameter_gentimeslot_mutant_active     <- rep(43, dim(parameter_table)[1] )
# duplet->triplet multiplier

parameter_table$parameter_AAAAAAAG_multiplier <- rep( 10^4, dim(parameter_table)[1] )

# mutiplier for powers definition
parameter_table$parameter_baseline_multiplier <- c(1,  seq(from = 100, to = 900, by  = 100) )
parameter_table$parameter_baseline_multiplier[length(parameter_table$parameter_baseline_multiplier)] <- 1000


# Overall mutation rate.
# Here the mutation rate per generation that will affect the overall mutation rate that is used annualy i.e. single time slot corresponds to 1 YEAR.
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000461
parameter_mutations_per_generation <- 5.2 * 10^{-10}

# Sugested numbe of base pairs in Antrax genome
# https://www.sciencedirect.com/science/article/pii/S151783821730206X
# taking 5,300,000 like Dr. Michael Norris has suggested.


# overall mutations rate
parameter_table$parameter01 <-  rep( parameter_mutations_per_generation , dim(parameter_table)[1] )



# Definging whether mutation interval is fixed or random.
parameter_table$parameter_interval_type <- rep( "fixed",  dim(parameter_table)[1] )

# defining number of time slots.
parameter_table$no_of_time_slots <- rep( 600,  dim(parameter_table)[1] )
                                             

# A
parameter_table$parameter02ac <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02ag <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02at <-  rep(0.25, dim(parameter_table)[1] )
# C
parameter_table$parameter02ca <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02cg <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02ct <-  rep(0.25, dim(parameter_table)[1] )
# G
parameter_table$parameter02ga <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02gc <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02gt <-  rep(0.25, dim(parameter_table)[1] )
# T
parameter_table$parameter02ta <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02tc <-  rep(0.25, dim(parameter_table)[1] )
parameter_table$parameter02tg <-  rep(0.25, dim(parameter_table)[1] )


# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000461
parameter_table$parameter03 <- rep( parameter_mutations_per_generation , dim(parameter_table)[1] )

# Picking wich one to insert probabilities.
parameter_table$parameter04a <- rep( 1/4, dim(parameter_table)[1] )
parameter_table$parameter04t <- rep( 1/4, dim(parameter_table)[1] )
parameter_table$parameter04c <- rep( 1/4, dim(parameter_table)[1] )

# Dropping probabilities
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000461
parameter_table$parameter05a <- rep( parameter_mutations_per_generation , dim(parameter_table)[1] )
parameter_table$parameter05t <- rep( parameter_mutations_per_generation , dim(parameter_table)[1] )
parameter_table$parameter05c <- rep( parameter_mutations_per_generation , dim(parameter_table)[1] )
parameter_table$parameter05g <- rep( parameter_mutations_per_generation , dim(parameter_table)[1] )


# Value of the fixation for wildtype and for mutant. It is assumed to be 1 for now.
# Active
parameter_table$parameter_pfix_mult_wildtype_active <- 1
parameter_table$parameter_pfix_mult_mutant_active   <- 1
# Sporadic
parameter_table$parameter_pfix_mult_wildtype_sporadic <- 1
parameter_table$parameter_pfix_mult_mutant_sporadic   <- 1


# Writing the results table into the file.
parameter_table_updated_relative_path <- "Design/parameters_table_extended_updated_final_fixed.csv"
write.csv(parameter_table, file = parameter_table_updated_relative_path, quote = FALSE, row.names = FALSE)

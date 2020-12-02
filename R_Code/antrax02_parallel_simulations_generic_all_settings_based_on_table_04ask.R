# // © Alexander Kirpich akirpich@gsu.edu
# The goal of this code is to write a mutations in Anthrax bacteria over time.
# Here the parallel scenario i is considered and scenarioID is determined from the design table.
rm(list=ls(all=TRUE))

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



# Reading parameters from the command line.
# This piece of code is used to aquire parameters from command line/in batch mode via sbatch submission of HiperGator.
args <- commandArgs(TRUE)
print(args)
current_scenario <- as.integer(args[1])
work_directory_path <- args[2]
parameter_table_relative_path <- args[3]
current_simulation_dataset_id <- as.integer(args[4])


# Debugging step to run on local machine instead instead of the code right above used for HiPer Gator.
# current_scenario <- 1
# current_scenario <- 10
# current_scenario <- 30
# work_directory_path  <- "E:/Alexander/UF Research/2019 - Cummings Norris Antrax"
# work_directory_path  <- "C:/Users/akirpich/Dropbox/GSU Research/2019 - Cummings Norris Antrax"
# parameter_table_relative_path <- "Design/parameters_table_extended_updated_final_fixed.csv"
# current_simulation_dataset_id <- 1
# current_simulation_dataset_id <- 595


# Setting up the working directory.
# This path is passed via work_directory_path variable that is passed via command line parameter on HiperGator.
setwd(work_directory_path)
# Extra check
# getwd()

# Setting seed for reproducibility.
# This seed is determinede based on the current_simulation_dataset_id which is different for every dataset that is simulated.
set.seed( current_simulation_dataset_id * 1000000 )


# Part 01
# Reading the table of settings inside R.
table_of_settings  <- read.table( parameter_table_relative_path,   sep = ",", header = TRUE )
# Debugging Step
# Printing parameters table on the screen.
print(table_of_settings)
cat(" current_scenario = ", current_scenario, "\n", sep ="")
cat(" scenario_id = ", as.character(table_of_settings$scenario_id[current_scenario]), "\n", sep ="")

# Fix 2019.05.19
# Reading fasta file rlative parth from the table.
fasta_file_with_data_relative_path <- as.character(table_of_settings$source_fasta_path[current_scenario])


# Part 02
# Reading original data into R that is used to generate the synthetic population.
fasta_file_data <- readDNAStringSet( filepath = fasta_file_with_data_relative_path, format = "fasta" )
# Debugging Step
# Printing parameters table on the screen.
print(fasta_file_data)

# Converting the readDNAStringSet to the sequence of vectors.
dataframe_file_data <- data.frame(fasta_file_data)
dim(dataframe_file_data)
# Extractring sample names.
sample_names <- rownames(dataframe_file_data)

# Parsing the dataset into the list of character sequences.
list_of_real_data_sequences <- lapply( seq_len(nrow(dataframe_file_data)), function(i) unlist( strsplit(x = dataframe_file_data[i,], split = NULL) ) )



# Part 03
# Generating synthetic population using the data from Part 03.
# Extracting the population size from the table and saving it into parameter_n.
parameter_n <- table_of_settings$parameter_n[current_scenario]
# parameter_n
# Extracting the number of the discrete time slots from the parameters table and saving it into no_of_time_slots.
no_of_time_slots <- table_of_settings$no_of_time_slots[current_scenario]

# Fix 2019.04.16
# This is an old code. The number of replications is not used anymore since each replications is taking its own thread.
# Extracting the number of replication
# no_of_replications <- table_of_settings$no_of_replications[current_scenario]

# n_dna_length
# Extracting the n_dna_length from the parameters table and saving it into n_dna_length.
n_dna_length <- table_of_settings$n_dna_length[current_scenario]

# parameter_AAAAAAAG_multiplier
# Extracting the parameter_AAAAAAAG_multiplier from the parameters table and saving it into parameter_AAAAAAAG_multiplier.
parameter_AAAAAAAG_multiplier <- table_of_settings$parameter_AAAAAAAG_multiplier[current_scenario]

# Fix 2019.09.11
# parameter_baseline_multiplier
# Extracting the parameter_baseline_multiplier from the parameters table and saving it into parameter_baseline_multiplier.
parameter_baseline_multiplier <- table_of_settings$parameter_baseline_multiplier[current_scenario]


# parameter_interval_type
# Extracting the parameter_interval_type from the parameters table and saving it into parameter_interval_type.
parameter_interval_type <- as.character(table_of_settings$parameter_interval_type[current_scenario])

# parameter_interval_time_mean
# Extracting the parameter_interval_time_mean from the parameters table and saving it into parameter_interval_time_mean.
parameter_interval_time_mean <- table_of_settings$parameter_interval_time_mean[current_scenario]

# parameter_interval_time_sd
# Extracting the parameter_interval_time_sd from the parameters table and saving it into parameter_interval_time_sd.
parameter_interval_time_sd <- table_of_settings$parameter_interval_time_sd[current_scenario]

# parameter_gentimeslot_wildtype_sporadic
# Extracting the parameter_gentimeslot_wildtype_sporadic from the parameters table and saving it into parameter_gentimeslot_wildtype_sporadic.
parameter_gentimeslot_wildtype_sporadic <- table_of_settings$parameter_gentimeslot_wildtype_sporadic[current_scenario]

# parameter_gentimeslot_wildtype_active
# Extracting the parameter_gentimeslot_wildtype_active from the parameters table and saving it into parameter_gentimeslot_wildtype_active.
parameter_gentimeslot_wildtype_active <- table_of_settings$parameter_gentimeslot_wildtype_active[current_scenario]

# parameter_gentimeslot_mutant_sporadic
# Extracting the parameter_gentimeslot_mutant_sporadic from the parameters table and saving it into parameter_gentimeslot_mutant_sporadic.
parameter_gentimeslot_mutant_sporadic <- table_of_settings$parameter_gentimeslot_mutant_sporadic[current_scenario]

# parameter_gentimeslot_mutant_active
# Extracting the parameter_gentimeslot_mutant_active from the parameters table and saving it into parameter_gentimeslot_mutant_active.
parameter_gentimeslot_mutant_active <- table_of_settings$parameter_gentimeslot_mutant_active[current_scenario]

# Part 04
# Extracting the parameters from the table that will be used
names(table_of_settings)

scenario_id   <- as.character(table_of_settings$scenario_id[current_scenario])

parameter01   <- table_of_settings$parameter01[current_scenario]

parameter02ac <- table_of_settings$parameter02ac[current_scenario]
parameter02at <- table_of_settings$parameter02at[current_scenario]
parameter02ag <- 1 - (parameter02ac + parameter02at)

parameter02ta <- table_of_settings$parameter02ta[current_scenario]
parameter02tc <- table_of_settings$parameter02tc[current_scenario]
parameter02tg <- 1 - (parameter02ta + parameter02tc)

parameter02ca <- table_of_settings$parameter02ca[current_scenario]
parameter02ct <- table_of_settings$parameter02ct[current_scenario]
parameter02cg <- 1 - (parameter02ca + parameter02ct)

parameter02ga <- table_of_settings$parameter02ga[current_scenario]
parameter02gt <- table_of_settings$parameter02gt[current_scenario]
parameter02gc <- 1 - (parameter02ga + parameter02gt)

parameter03   <- table_of_settings$parameter03[current_scenario]

parameter04a  <- table_of_settings$parameter04a[current_scenario]
parameter04t  <- table_of_settings$parameter04t[current_scenario]    
parameter04c  <- table_of_settings$parameter04c[current_scenario]
parameter04g  <- 1 - (parameter04a + parameter04t + parameter04c)


parameter05a  <- table_of_settings$parameter05a[current_scenario] 
parameter05t  <- table_of_settings$parameter05t[current_scenario]
parameter05c  <- table_of_settings$parameter05c[current_scenario]
parameter05g  <- table_of_settings$parameter05g[current_scenario]

# Fix 2019.09.10
# Extracting fixation parameters from the table.
# Active
parameter_pfix_mult_wildtype_active   <- table_of_settings$parameter_pfix_mult_wildtype_active[current_scenario]
parameter_pfix_mult_mutant_active     <- table_of_settings$parameter_pfix_mult_mutant_active[current_scenario]
# Sporadic
parameter_pfix_mult_wildtype_sporadic <- table_of_settings$parameter_pfix_mult_wildtype_sporadic[current_scenario]
parameter_pfix_mult_mutant_sporadic   <- table_of_settings$parameter_pfix_mult_mutant_sporadic[current_scenario]




# Fix 2019.04.16
# Here we are using single thread for the single simulation
# The old piece of code is used below where current_replication is assinged tuccent job id 
current_replication <- current_simulation_dataset_id

# Debuggins step for individuals
# current_replication <-  1 


# Generating synthetic population indexes.
list_of_sample_indexes_of_sampled_population <- sample( x = c(1:dim(dataframe_file_data)[1]), size = parameter_n, replace = TRUE )
# Summaries of the just generated population
# unique(list_of_sample_indexes_of_sampled_population)
# length( unique(list_of_sample_indexes_of_sampled_population)  )
# table(unique(list_of_sample_indexes_of_sampled_population))

# Generating the synthetic population
list_of_synthetic_population_sequences <- list_of_real_data_sequences[list_of_sample_indexes_of_sampled_population]
dim( summary(list_of_synthetic_population_sequences) )
# list_of_synthetic_population_sequences
# Length of the sequence for a population element.
initial_sequence_length <- length( list_of_synthetic_population_sequences[[1]] )

# Generating synthetic population names
sample_names_synthetic    <- paste(sample_names[list_of_sample_indexes_of_sampled_population])
sample_numbers_synthetic <- paste( "sample_element_", seq(1:length(list_of_sample_indexes_of_sampled_population)), "_from_", sep =""  )
sample_names_synthetic_population <- paste(sample_numbers_synthetic, sample_names_synthetic, sep ="")


dim( summary(list_of_synthetic_population_sequences) )
# list_of_synthetic_population_sequences


# Generating the table that will store stop codon presence  over time.
# This table has no_of_time_slots columns for each individual unit in the synthetic population.
stop_codon_presence_over_time <- matrix( "X", nrow = parameter_n, ncol = no_of_time_slots )
dim(stop_codon_presence_over_time)
stop_codon_presence_over_time[ c(1:10), c(1:10) ]
# Postulating stop codons that will be used for matching.
codon_tag <- "TAG" 
codon_taa <- "TAA" 
codon_tga <- "TGA"
# Postulating coding codons starting points.
# Here the asssumtion is that the lenght of the original DNA piece can increase 5 fold so this is accounted for in the starting condon points generation.
codon_starting_points <- seq(from = 1, to = (initial_sequence_length*15), by = 3)
# codon_starting_points[1:10]
length(codon_starting_points)


# Fix 2019.09.09.
# Generating the table that will store wildtype/Mutant status over time.
# This table has no_of_time_slots columns for each individual unit in the synthetic population.
wildtype_mutant_status_over_time <- matrix( "X", nrow = parameter_n, ncol = no_of_time_slots )
dim(wildtype_mutant_status_over_time)
wildtype_mutant_status_over_time[ c(1:10), c(1:10) ]
# Postulating stop codons that will be used for matching.
codon_doubleAAAAAAAG <- "AAAAAAAGAAAAAAAG" 
codon_tripleAAAAAAAG <- "AAAAAAAGAAAAAAAGAAAAAAAG" 
# Defining length for both double and triplet
# Double
length_codon_doubleAAAAAAAG <- length( unlist(strsplit(x = codon_doubleAAAAAAAG, split = NULL)))
# Triple
length_codon_tripleAAAAAAAG <- length( unlist(strsplit(x = codon_tripleAAAAAAAG, split = NULL)))




# Fix 2019.05.04
# Generating the table that will store hamming distance over time.
# This table has no_of_time_slots columns for each individual unit in the synthetic population.
hamming_distance_over_time <- matrix( NA, nrow = parameter_n, ncol = no_of_time_slots )
dim(hamming_distance_over_time)
hamming_distance_over_time[ c(1:10), c(1:10) ]
sum( is.na(hamming_distance_over_time) )


# Fix 2019.09.11.
# This part can be discussed with Michael Norris in more details!!!
# It can be adjusted in the future if necessary.

# Generating fixation parameters for wildtype/mutant and active/sporadic states.
# Sporadic
parameter_pfix_wildtype_sporadic <-  parameter_pfix_mult_wildtype_sporadic
parameter_pfix_mutant_sporadic   <-  parameter_pfix_mult_mutant_sporadic
# Active
parameter_pfix_wildtype_active   <-  parameter_pfix_mult_wildtype_active
parameter_pfix_mutant_active     <-  parameter_pfix_mult_mutant_active




# Part 05 looping over all individuals over time
# This piece of code loops individual by individual in a synthetic population and over all time slots.

# Duplicating a copy so that we can track the original version
list_of_synthetic_population_sequences_original <- list_of_synthetic_population_sequences



# Fix 2019.09.09
# Here we fill those year that have active status.
# This filling depends on whether the filling mechanism is (fixed) or (random).
# Regardless of the algorithm in the end list_of_active_years variable contains the 
# (fixed)
if ( parameter_interval_type == "fixed" ) 
{ 
  # Generating list of active years based on the fixed assumption
  list_of_active_years <- seq( from = 1, to = no_of_time_slots, by = parameter_interval_time_mean )
}  
# (random)
if ( parameter_interval_type == "random" ) 
{
  # Generating list of active periods b/w years based on the mean and sd
  list_of_normal_rv_outcomes <- rnorm( n = (no_of_time_slots/parameter_interval_time_mean) * 100 , mean = parameter_interval_time_mean, sd = parameter_interval_time_sd )
  # Rounding this list since integer increments are of interest
  list_of_normal_rv_outcomes_rounded <- round( list_of_normal_rv_outcomes )
  # Getting indexes using cummulative sum funtion
  list_of_normal_rv_outcomes_rounded_cumsums <- cumsum( x = list_of_normal_rv_outcomes_rounded )
  # Keeping only those indexes which are below no_of_time_slots variable.
  list_of_active_years <- list_of_normal_rv_outcomes_rounded_cumsums[ list_of_normal_rv_outcomes_rounded_cumsums <= no_of_time_slots ]
}  




# Looping over individuals in the current simulation
for( i in c(1:parameter_n) )
{
  
  # Debuggins step for individuals
  # i <-  1 

  # Determinining original sequence before manipulations.
  original_sequence_length_i <- length(unlist(list_of_synthetic_population_sequences_original[[i]]))
  
  
  # Extracting current indiviudual in a synthetic population that will be modified in process.
  current_synthetic_sequence <- list_of_synthetic_population_sequences_original[[i]]

  # Fix 2019.05.04
  # Saving the original sequence. It will be comared with the constantly evolving current_synthetic_sequence from previous step
  original_synthetic_sequence_string <- paste( list_of_synthetic_population_sequences_original[[i]], collapse = '' )
  original_synthetic_sequence_string_dna <- DNAString(original_synthetic_sequence_string)

  
  # Debuggin step
  # no_of_time_slots <- 3
  # no_of_time_slots <- 50
  
  # Looping over all tieme slot for a given individual i.
  for( current_time_slot in c(1:no_of_time_slots) )
  {
    # Debugging step
    # current_time_slot <- 1
    # current_time_slot <- 2
    
    
    # Fix 2019.09.09
    # Determining stop codon presence  at the current_time_slot and filling the result into stop_codon_presence_over_time  
    
    # Converting it to a string.
    current_synthetic_sequence_string <- paste( current_synthetic_sequence, collapse = '' )


    # Matching codons with the string for all the stopping codons.
    codon_tag_matches <- gregexpr( pattern = codon_tag, text = current_synthetic_sequence_string )
    codon_taa_matches <- gregexpr( pattern = codon_taa, text = current_synthetic_sequence_string )
    codon_tga_matches <- gregexpr( pattern = codon_tga, text = current_synthetic_sequence_string )
    # Merging starting points for the coding codons.
    codon_matches_current_synthetic_sequence <- c( as.vector(codon_tag_matches[[1]]),  as.vector(codon_taa_matches[[1]]), as.vector(codon_tga_matches[[1]]) ) 
    # Checking whether the codon starting points are in in the codon_starting_points sequence.
    # If it is there the current sequence is coded as Inaba otherwise it is coded as Ogawa.
    sum_codon_matches_current_synthetic_sequence <- sum( codon_matches_current_synthetic_sequence %in% codon_starting_points )
    # Doing the appropriate assingment.
    if ( sum_codon_matches_current_synthetic_sequence > 0  )  { stop_codon_presence_over_time[i, current_time_slot] <- "P" }
    if ( sum_codon_matches_current_synthetic_sequence == 0 )  { stop_codon_presence_over_time[i, current_time_slot] <- "A" }
    
    
    # Fix 2019.09.09.
    # Matching codons with the string forthe stopping codons.
    # Double
    codon_doubleAAAAAAAG_matches <- gregexpr( pattern = codon_doubleAAAAAAAG, text = current_synthetic_sequence_string )
    # Triple
    codon_tripleAAAAAAAG_matches <- gregexpr( pattern = codon_tripleAAAAAAAG, text = current_synthetic_sequence_string )
    
    # Merging starting points for the coding codons.
    codon_matches_current_synthetic_sequence <- c( as.vector(codon_tripleAAAAAAAG_matches[[1]])  ) 
    # Checking whether the codon starting points are in in the codon_starting_points sequence.
    # If it is there the current sequence is coded as Mutant otherwise it is coded as wildtype.
    sum_codon_matches_current_synthetic_sequence <- sum( codon_matches_current_synthetic_sequence %in% codon_starting_points )
    # Doing the appropriate assingment.
    if ( sum_codon_matches_current_synthetic_sequence > 0  )  { wildtype_mutant_status_over_time[i, current_time_slot] <- "M" }
    if ( sum_codon_matches_current_synthetic_sequence == 0 )  { wildtype_mutant_status_over_time[i, current_time_slot] <- "W" }
    
    
    
    # Fix 2019.09.09.
    # Computing Hamming distance b/w the original sequence and the current sequence.
    # Generating slots to save hamming distance. 
    list_of_time_slos_to_save <- c(1,c(1:floor(no_of_time_slots/25) * 25))
    if ( current_time_slot %in% list_of_time_slos_to_save  )
    {
      
      # At first alignmnet has to be performed.
      # Generating DNA string type object.
      current_synthetic_sequence_string_dna <- DNAString(current_synthetic_sequence_string)
      # Performing alignment current_synthetic_sequence_string_dna to the original_synthetic_sequence_string_dna
      current_alignment_output <- pairwiseAlignment(subject = original_synthetic_sequence_string_dna , pattern = current_synthetic_sequence_string_dna, type = "global-local",
                                                    gapOpening = 0, gapExtension = 0 )
      # Extracting the results of the alignment.
      original_synthetic_sequence_string_dna_aligned <- as.character( attr(current_alignment_output, "subject") )
      current_synthetic_sequence_string_dna_aligned  <- as.character( attr(current_alignment_output, "pattern") )
      # Converting strings back to the array of characters.
      original_synthetic_sequence_dna_aligned <- unlist( strsplit(x = original_synthetic_sequence_string_dna_aligned, split = NULL) )
      current_synthetic_sequence_dna_aligned  <- unlist( strsplit(x = current_synthetic_sequence_string_dna_aligned,  split = NULL) )
      # computing current Hamming distance and saving it into the matrix.
      hamming_distance_over_time[i, current_time_slot] <- original_sequence_length_i - sum( current_synthetic_sequence_dna_aligned == original_synthetic_sequence_dna_aligned )      
    }  

    
    
    
    # Performing 7 d) from the draft.
    # Here we decide whether to keep the mutations that will happen.
    # The probability of mutations depends on whether it is
    # -- Widtype or mutant at the CURRENT time slot.
    # -- Whether it is sporadic year or active year from the 
    # fixation_probability is defined below.
    # Wildtype
    if ( ( wildtype_mutant_status_over_time[i, current_time_slot] == "W" ) && ( !(current_time_slot %in% list_of_active_years) ) )  
    {

      # Fixation probability
      fixation_probability <-  parameter_pfix_wildtype_sporadic 
      # Mutation rate multiplier
      parameter_current_probability_multiplier <- parameter_gentimeslot_wildtype_sporadic
    }
    if ( ( wildtype_mutant_status_over_time[i, current_time_slot] == "W" ) && (   current_time_slot %in% list_of_active_years  ) )  
    { 
      
      # Fixation probability
      fixation_probability <-  parameter_pfix_wildtype_active  
      # Mutation rate multiplier
      parameter_current_probability_multiplier <- parameter_gentimeslot_wildtype_active
    }
    # Mutant
    if ( ( wildtype_mutant_status_over_time[i, current_time_slot] == "M" ) && ( !(current_time_slot %in% list_of_active_years) ) )  
    { 
      
      # Fixation probability
      fixation_probability <-  parameter_pfix_mutant_sporadic 
      # Mutation rate multiplier
      parameter_current_probability_multiplier <- parameter_gentimeslot_mutant_sporadic
    }
    if ( ( wildtype_mutant_status_over_time[i, current_time_slot] == "M" ) && (   current_time_slot %in% list_of_active_years  ) )  
    { 
      
      # Fixation probability
      fixation_probability <-  parameter_pfix_mutant_active  
      # Mutation rate multiplier
      parameter_current_probability_multiplier <- parameter_gentimeslot_mutant_active
    }
    

    
    # Determining whether fixation is performed or not.
    fixation_decision  <-  rbinom( n = 1, size = 1, prob = fixation_probability )
    
    # Debugging step
    # fixation_decision <- 1
    
    
    # Performing steps 5 ), 5 b) and 5 c) only if fixation_decision == 1 i.e. if there is a decision to fix the mutation that is happening at the given time step.
    # If there is no decision to fix the mutation then the simulation is going successfuly to the next step without changes.
    if ( fixation_decision == 1 )
    { 
      
      
      # Fix 2019.09.09.
      # Mutations are performed for duplet -> triplet regardless of the subsequent steps.
      # This step is performed only for wildtype Anthrax, therefore it is checked whether it is.
      # We are mutating the wildtype to mutant if it is not mutant already.
      # In this case we are trying to move it from Wildtype to mutant with probability equal to parameter_AAAAAAAG_multiplier * <probabily of insertion>.
      # This is performed only if the current dna in current time is still wildtype and NOT mutant.
      
      # In the begiing we are duplicating current_synthetic_sequence into current_synthetic_sequence_afterAAAAAAAG.
      # This current_synthetic_sequence_afterAAAAAAAG will be updated only if there are conditions below o/w it will stay current synthetic sequence.
      current_synthetic_sequence_afterAAAAAAAG <- current_synthetic_sequence

      if ( wildtype_mutant_status_over_time[i, current_time_slot] == "W" ) 
      {
        
        current_mutation_AAAAAAAG_status <- rbinom( n = 1, size = 1, prob = parameter03 * parameter_baseline_multiplier * parameter_AAAAAAAG_multiplier )
        
        # Debugign step
        # current_mutation_AAAAAAAG_status <- 1
        
        # If the decision is to move from wildtype to mutant this is performed.
        if ( current_mutation_AAAAAAAG_status == 1 ) 
        {
          
          starting_index_doubleAAAAAAAG <- gregexpr( pattern = codon_doubleAAAAAAAG, text = paste( current_synthetic_sequence, collapse = '' ) )[[1]][1]
          
          indexes_doubleAAAAAAAG <- c(starting_index_doubleAAAAAAAG:(starting_index_doubleAAAAAAAG + length_codon_doubleAAAAAAAG - 1 ))
          indexes_singleAAAAAAAG <- indexes_doubleAAAAAAAG[ (1:(length(indexes_doubleAAAAAAAG)/2)) ]
          
          # Extracting sequence to insert
          sequence_to_insert <- current_synthetic_sequence[ indexes_singleAAAAAAAG ]
          
          # Parsing the original sequnce 
          current_synthetic_sequence_parse_before  <- current_synthetic_sequence[c(1:(starting_index_doubleAAAAAAAG-1))]
          current_synthetic_sequence_parse_after   <- current_synthetic_sequence[c(starting_index_doubleAAAAAAAG:length(current_synthetic_sequence))]
          # Recombining the sequence
          current_synthetic_sequence_afterAAAAAAAG <- c(current_synthetic_sequence_parse_before, sequence_to_insert, current_synthetic_sequence_parse_after)
          
          # Debugging step
          # current_synthetic_sequence_mutation[ c(starting_index_doubleAAAAAAAG:(starting_index_doubleAAAAAAAG+30)) ]
          
        }
        
      }  
      

      
      # Performing 5 a) from the draft.
      
      # Fix 2019.09.09
      # Excluding duplets and triplets from further mutations
      
      # Searching for patterns
      # Double
      codon_doubleAAAAAAAG_matches_current <- gregexpr( pattern = codon_doubleAAAAAAAG, text = paste( current_synthetic_sequence_afterAAAAAAAG, collapse = '' )  )
      # Triple
      codon_tripleAAAAAAAG_matches_current <- gregexpr( pattern = codon_tripleAAAAAAAG, text = paste( current_synthetic_sequence_afterAAAAAAAG, collapse = '' )  )

      # Indentifying indexes
      # Double
      codon_doubleAAAAAAAG_matches_current_indexes <- c( as.vector(codon_doubleAAAAAAAG_matches_current[[1]]):(as.vector(codon_doubleAAAAAAAG_matches_current[[1]])+length_codon_doubleAAAAAAAG-1) ) 
      # Triple (if there is Triple)
      if (as.vector(codon_tripleAAAAAAAG_matches_current[[1]]) > 0 )
      {
        codon_tripleAAAAAAAG_matches_current_indexes <- c( as.vector(codon_tripleAAAAAAAG_matches_current[[1]]):(as.vector(codon_tripleAAAAAAAG_matches_current[[1]])+length_codon_tripleAAAAAAAG-1) ) 
      }  else { codon_tripleAAAAAAAG_matches_current_indexes <- NULL  }  
      
      # Generating combined indexes
      codon_AAAAAAAG_matches_current_indexes <- union(codon_doubleAAAAAAAG_matches_current_indexes, codon_tripleAAAAAAAG_matches_current_indexes)
      
      
      # Backing up current sequence before it is mutated.
      current_synthetic_sequence_mutation <- current_synthetic_sequence_afterAAAAAAAG

      
      # Generating mutations for each character in the sequence.
      current_mutated_characters <- rbinom( n = length(current_synthetic_sequence_afterAAAAAAAG), size = 1, prob = parameter01 * parameter_baseline_multiplier * parameter_current_probability_multiplier )
      # Fix 2019.09.10
      # Dropping those that belong to duplet/triplet codon codon_AAAAAAAG. They will not be mutated regardless.
      current_mutated_characters_places <- setdiff( which( current_mutated_characters == 1 ) , codon_AAAAAAAG_matches_current_indexes )
      # Debugging
      # length(current_mutated_characters_places)
            
      
      # Determining the actual mutation shift.  
      for( current_mutation_index in current_mutated_characters_places )
      {
        
        # Debugging step
        # current_mutation_index <- current_mutated_characters_places[1]
        
        # Determining mutationf for each A, C, T, and G.
        # Working with A  
        if ( current_synthetic_sequence_mutation[current_mutation_index] == "A" )
        {
          # Generating the current mutation for A.   
          current_mutation_outcome  <-  rmultinom( n = 1,  size = 1,  prob = c(parameter02ac, parameter02at, parameter02ag)  )
          
          # Performing the corresponding change from A.
          if ( current_mutation_outcome[1] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "C" }
          if ( current_mutation_outcome[2] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "T" }
          if ( current_mutation_outcome[3] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "G" }
          
          # End of -> if ( current_synthetic_sequence_mutation[current_mutation_index] == "A" )  
        }
        
        
        # Working with T  
        if ( current_synthetic_sequence_mutation[current_mutation_index] == "T" )
        {
          # Generating the current mutation for T.   
          current_mutation_outcome  <-  rmultinom( n = 1,  size = 1,  prob = c(parameter02ta, parameter02tc, parameter02tg)  )
          
          # Performing the corresponding change from T.
          if ( current_mutation_outcome[1] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "A" }
          if ( current_mutation_outcome[2] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "C" }
          if ( current_mutation_outcome[3] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "G" }
          
          # End of -> if ( current_synthetic_sequence_mutation[current_mutation_index] == "T" )  
        }
        
        
        # Working with C  
        if ( current_synthetic_sequence_mutation[current_mutation_index] == "C" )
        {
          # Generating the current mutation for C.   
          current_mutation_outcome  <-  rmultinom( n = 1,  size = 1,  prob = c(parameter02ca, parameter02ct, parameter02cg)  )
          
          # Performing the corresponding change from C.
          if ( current_mutation_outcome[1] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "A" }
          if ( current_mutation_outcome[2] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "T" }
          if ( current_mutation_outcome[3] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "G" }
          
          # End of -> if ( current_synthetic_sequence_mutation[current_mutation_index] == "C" )  
        }
        
        
        # Working with G  
        if ( current_synthetic_sequence_mutation[current_mutation_index] == "G" )
        {
          # Generating the current mutation for G.   
          current_mutation_outcome  <-  rmultinom( n = 1,  size = 1,  prob = c(parameter02ga, parameter02gc, parameter02gt)  )
          
          # Performing the corresponding change from C.
          if ( current_mutation_outcome[1] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "A" }
          if ( current_mutation_outcome[2] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "C" }
          if ( current_mutation_outcome[3] == 1 )  {  current_synthetic_sequence_mutation[current_mutation_index] <- "T" }
          
          # End of -> if ( current_synthetic_sequence_mutation[current_mutation_index] == "G" )  
        }
        
        # End of -> for( current_mutation_index in current_mutated_characters_places )  
      }  
      
      
      
      
      # Performing 5 b) from the draft.
      # Determining the current length
      
      # Generating insertion slots with the numbers length(current_synthetic_sequence_mutation) + 1 since the number of spots is larger by 1 from length(current_synthetic_sequence_mutation)
      # that was inhereted from the previous step that simulated mutations.
      current_inserted_characters <- rbinom( n = length(current_synthetic_sequence_mutation) + 1, size = 1, prob = parameter03 * parameter_baseline_multiplier * parameter_current_probability_multiplier )
      # Debugging step
      # current_inserted_characters <- c( rep(1, (length(current_synthetic_sequence_mutation) + 1) * 2/3 ), rep(0, (length(current_synthetic_sequence_mutation) + 1) * 1/3 ) )
      
      # Fix 2009.09.10
      # Dropping those that belong to duplet/triplet codon codon_AAAAAAAG. They will not undergo insertations regardless.
      # The only ones that cannot be inserted are the ones that are WITHIN the codon, therefore codon_AAAAAAAG_matches_current_indexes has to be adjusted.
      # + 1 has to be added to the current index to exclude before indexes.
      current_inserted_characters_places <- setdiff( which( current_inserted_characters == 1 ) , codon_AAAAAAAG_matches_current_indexes[-1] )
      # Debugging
      # length(current_mutated_characters_places)

      
      # The imputations have to be performed from the last index since the end of the array will change on every iteration.
      # sort(current_inserted_characters_places, decreasing = TRUE)
      
      
      # Backing up current sequence before it is modified by insertions.
      current_synthetic_sequence_insertion <- current_synthetic_sequence_mutation
      
      # Determining the actual insertion characters.  
      for( current_insertion_index in sort(current_inserted_characters_places, decreasing = TRUE) )
      {
        
        # Debugging step
        # current_insertion_index <- sort(current_inserted_characters_places, decreasing = TRUE)[1]
        # current_insertion_index <- sort(current_inserted_characters_places, decreasing = TRUE)[2]      
        
        # Considering border cases separately since thye are inserted differently
        # First index
        if ( current_insertion_index == 1 )
        {
          # Parsing the array to be able to perform the insertation inbetween the pieces at  current_insertion_index.
          current_synthetic_sequence_part01 <- NULL
          current_synthetic_sequence_part02 <- current_synthetic_sequence_insertion
        }
        # Last index
        if ( current_insertion_index == (length(current_synthetic_sequence_insertion)+1) )
        {
          # Parsing the array to be able to perform the insertation inbetween the pieces at  current_insertion_index.
          current_synthetic_sequence_part01 <- current_synthetic_sequence_insertion
          current_synthetic_sequence_part02 <- NULL
        }
        
        # Considering middle cases. 
        if ( (current_insertion_index>1) && (current_insertion_index<length(current_synthetic_sequence_insertion)+1) )
        {
          # Parsing the array to be able to perform the insertation inbetween the pieces at  current_insertion_index.
          current_synthetic_sequence_part01 <- current_synthetic_sequence_insertion[ c(1:(current_insertion_index-1)) ]
          current_synthetic_sequence_part02 <- current_synthetic_sequence_insertion[ c(current_insertion_index:length(current_synthetic_sequence_insertion)) ]
        }
        
        
        # Debugging piece
        # cat("current_insertion_index ->", current_insertion_index, "\n")
        # cat("01 ->", length(current_synthetic_sequence_part01), "\n")
        # cat("02 ->", length(current_synthetic_sequence_part02), "\n")
        # cat("01+02 ->", length(current_synthetic_sequence_part01) + length(current_synthetic_sequence_part02), "\n")
        # cat("BEFORE: length( current_synthetic_sequence_insertion) ->",  length(current_synthetic_sequence_insertion), "\n")
        
        # Determining insertion character for each insertion.
        
        # Generating the current insertion character.   
        current_insertion_outcome_character <- rmultinom( n = 1,  size = 1,  prob = c(parameter04a, parameter04t, parameter04c, parameter04g)  )
        
        # Performing the corresponding insertation.
        if ( current_insertion_outcome_character[1] == 1 )  {  current_synthetic_sequence_insertion <-  c(current_synthetic_sequence_part01, "A", current_synthetic_sequence_part02) }
        if ( current_insertion_outcome_character[2] == 1 )  {  current_synthetic_sequence_insertion <-  c(current_synthetic_sequence_part01, "T", current_synthetic_sequence_part02) }
        if ( current_insertion_outcome_character[3] == 1 )  {  current_synthetic_sequence_insertion <-  c(current_synthetic_sequence_part01, "C", current_synthetic_sequence_part02) }
        if ( current_insertion_outcome_character[4] == 1 )  {  current_synthetic_sequence_insertion <-  c(current_synthetic_sequence_part01, "G", current_synthetic_sequence_part02) }
        
        
        # Debugging piece
        # cat("AFTER: length( current_synthetic_sequence_insertion) ->",  length(current_synthetic_sequence_insertion), "\n\n")
        
        # End of -> for( current_insertion_index in current_inserted_characters_places )
      }  
      
      
      
      
      
      
      # Performing 5 c) from the draft
      
      # Fix 2019.09.03
      # Removng empty elements that appeared on the previous step.
      current_synthetic_sequence_insertion <- current_synthetic_sequence_insertion[current_synthetic_sequence_insertion != ""]
      
      # Determining the current length
      current_sequence_length <- length(current_synthetic_sequence_insertion)
      
      
      # Generating deletion slots for current_sequence_length numbers.
      # Deletion probabilities are different for all characters i.e. they are differetn for A, T, C and G.
      # A 
      which_character_a <-  which(current_synthetic_sequence_insertion == "A")
      current_deleted_characters_a <- rbinom( n = length(which_character_a), size = 1, prob = parameter05a * parameter_baseline_multiplier * parameter_current_probability_multiplier )
      current_deleted_characters_a_local_places  <- which( current_deleted_characters_a == 1 )
      current_deleted_characters_a_global_places <- which_character_a[current_deleted_characters_a_local_places]
      # T 
      which_character_t <-  which(current_synthetic_sequence_insertion == "T")
      current_deleted_characters_t <- rbinom( n = length(which_character_t), size = 1, prob = parameter05t * parameter_baseline_multiplier * parameter_current_probability_multiplier )
      current_deleted_characters_t_local_places  <- which( current_deleted_characters_t == 1 )
      current_deleted_characters_t_global_places <- which_character_t[current_deleted_characters_t_local_places]
      # C 
      which_character_c <-  which(current_synthetic_sequence_insertion == "C")
      current_deleted_characters_c <- rbinom( n = length(which_character_c), size = 1, prob = parameter05c * parameter_baseline_multiplier * parameter_current_probability_multiplier )
      current_deleted_characters_c_local_places  <- which( current_deleted_characters_c == 1 )
      current_deleted_characters_c_global_places <- which_character_c[current_deleted_characters_c_local_places]
      # G 
      which_character_g <-  which(current_synthetic_sequence_insertion == "G")
      current_deleted_characters_g <- rbinom( n = length(which_character_g), size = 1, prob = parameter05g * parameter_baseline_multiplier * parameter_current_probability_multiplier )
      current_deleted_characters_g_local_places  <- which( current_deleted_characters_g == 1 )
      current_deleted_characters_g_global_places <- which_character_g[current_deleted_characters_g_local_places]
      
      # Combining characters to delete.
      current_deleted_characters_all_global_places_preliminary_with_codonAAAAAAAG <- c( current_deleted_characters_a_global_places, current_deleted_characters_t_global_places,
                                                                                        current_deleted_characters_c_global_places, current_deleted_characters_g_global_places )
      
      
      
      # Fix 2019.09.10
      # Excluding duplets and triplets from further deletion
      
      # Searching for patterns
      # Double
      codon_doubleAAAAAAAG_matches_current_deletion <- gregexpr( pattern = codon_doubleAAAAAAAG, text = paste( current_synthetic_sequence_insertion, collapse = '' )  )
      # Triple
      codon_tripleAAAAAAAG_matches_current_deletion <- gregexpr( pattern = codon_tripleAAAAAAAG, text = paste( current_synthetic_sequence_insertion, collapse = '' )  )
      
      # Indentifying indexes
      # Double
      codon_doubleAAAAAAAG_matches_current_indexes_deletion <- c( as.vector(codon_doubleAAAAAAAG_matches_current_deletion[[1]]):(as.vector(codon_doubleAAAAAAAG_matches_current_deletion[[1]])+length_codon_doubleAAAAAAAG-1) ) 
      # Triple (if there is Triple)
      if (as.vector(codon_tripleAAAAAAAG_matches_current_deletion[[1]]) > 0 )
      {
        codon_tripleAAAAAAAG_matches_current_indexes_deletion <- c( as.vector(codon_tripleAAAAAAAG_matches_current_deletion[[1]]):(as.vector(codon_tripleAAAAAAAG_matches_current_deletion[[1]])+length_codon_tripleAAAAAAAG-1) ) 
      }  else { codon_tripleAAAAAAAG_matches_current_indexes_deletion <- NULL  }  
      
      # Generating combined indexes
      codon_AAAAAAAG_matches_current_indexes_deletion <- union(codon_doubleAAAAAAAG_matches_current_indexes_deletion, codon_tripleAAAAAAAG_matches_current_indexes_deletion)

      # Finally combining characters to delete and excluding the 
      # Combining characters to delete.
      current_deleted_characters_all_global_places <-  setdiff( current_deleted_characters_all_global_places_preliminary_with_codonAAAAAAAG, codon_AAAAAAAG_matches_current_indexes_deletion  )
      

      
      # Fix 2019.08.30.
      # Identifying characters to keep. 
      if ( length(current_synthetic_sequence_insertion) == 0 )
      {  
        current_kept_characters_all_global_places    <- NULL
      }
      if ( length(current_synthetic_sequence_insertion) > 0 )
      {  
        current_kept_characters_all_global_places    <- setdiff( c(1:length(current_synthetic_sequence_insertion)), current_deleted_characters_all_global_places )
      }
      
      # Extra check
      # sum( ! sort( union( current_deleted_characters_all_global_places, current_kept_characters_all_global_places) ) == c(1:length(current_synthetic_sequence_insertion)) )
      
      # Fix 2019.08.30.
      # Identifying characters to keep. 
      if ( length(current_synthetic_sequence_insertion) == 0 )
      {  
        current_synthetic_sequence_deletion  <- ""
      }
      if ( length(current_synthetic_sequence_insertion) > 0 )
      {  
        # Backing up current sequence before it is modified by deletion.
        current_synthetic_sequence_deletion <- unlist(strsplit(x = current_synthetic_sequence_insertion, split = NULL)[current_kept_characters_all_global_places])
        
        # Fix 2019.09.03
        # Checking if the resulting sequecne is NULL and has lenght zero. In this case 
        if ( length(current_synthetic_sequence_deletion) == 0 )
        {  
          current_synthetic_sequence_deletion <- ""
        }  
      }
      
      
      # Saving the current synthetic sequence that will be passed to the next time slot.
      current_synthetic_sequence <- current_synthetic_sequence_deletion
      
      # End of -> if ( fixation_decision == 1 )  
    }    
    
    
    
    # End of -> for( current_time_slot in c(1:no_of_time_slots) )  
  }  
  
  # Returning the result to the original list.
  # Here the original sewuence is overwrittent.
  list_of_synthetic_population_sequences[[i]] <- current_synthetic_sequence
  
  # Fix 2019.08.29.
  # Debuggin step
  cat("Individual Sequence -> ", i, " Simulated Fine! \n")
  
  
  # End of -> for( i in c(1:parameter_n) )  
}  


# Generating a version that can be exported into the fasta file for futher reading.
list_of_synthetic_population_sequences_strings <- list_of_synthetic_population_sequences
# Same for the original data before changes
list_of_synthetic_population_sequences_original_strings <- list_of_synthetic_population_sequences_original


# overriding each sequence to create the strings version. 
for( i in c(1:parameter_n) )
{
  
  # Debuggins step for individuals
  # i <-  1 
  
  list_of_synthetic_population_sequences_strings[[i]] <- paste(list_of_synthetic_population_sequences[[i]], collapse = '')
  
  # Same for the original data before changes.
  list_of_synthetic_population_sequences_original_strings[[i]] <- paste(list_of_synthetic_population_sequences_original[[i]], collapse = '')
  
  # End of -> for( i in c(1:parameter_n) )  
}  




# Writing out the results of the simulations
# TO do that each file/object has to inlcude the run sequence id which is saved in the curren_replication


# ---- Creating object with the correct namings ----

# Generating object for list of arrays
eval( parse( text=paste( "list_of_synthetic_population_sequences_scenario_id_", scenario_id, "_replication_", current_replication, " <- list_of_synthetic_population_sequences", sep ="") )  )
# Generating object for list of strings
eval( parse( text=paste( "list_of_synthetic_population_sequences_strings_scenario_id_", scenario_id, "_replication_", current_replication, " <- list_of_synthetic_population_sequences_strings", sep ="") )  )

# Same for the original array
# Generating object for list of arrays
eval( parse( text=paste( "list_of_synthetic_population_sequences_original_scenario_id_", scenario_id, "_replication_", current_replication, " <- list_of_synthetic_population_sequences_original", sep ="") )  )
# Generating object for list of strings
eval( parse( text=paste( "list_of_synthetic_population_sequences_original_strings_scenario_id_", scenario_id, "_replication_", current_replication, " <- list_of_synthetic_population_sequences_original_strings", sep ="") )  )


# ---- Creating path with the correct namings ----

# Generating saving path for list of arrays
list_of_synthetic_population_sequences_current_replication_path <- paste("R_Data/list_of_synthetic_population_sequences_scenario_id_", scenario_id, "_replication_", current_replication, ".RData", sep ="")
# Generating saving path for list of strings
list_of_synthetic_population_sequences_strings_current_replication_path <- paste("R_Data/list_of_synthetic_population_sequences_strings_scenario_id_", scenario_id, "_replication_", current_replication, ".RData", sep ="")

# Same for the original array
# Generating saving path for list of arrays
list_of_synthetic_population_sequences_original_current_replication_path <- paste("R_Data/list_of_synthetic_population_sequences_original_scenario_id_", scenario_id, "_replication_", current_replication, ".RData", sep ="")
# Generating saving path for list of strings
list_of_synthetic_population_sequences_original_strings_current_replication_path <- paste("R_Data/list_of_synthetic_population_sequences_original_strings_scenario_id_", scenario_id, "_replication_", current_replication, ".RData", sep ="")


# ---- Creating path for the fasta and summary tables ----

# Generating saving path for fasta file
list_of_synthetic_population_sequences_fasta_current_replication_path   <- paste("R_Output/list_of_synthetic_population_sequences_fasta_scenario_id_", scenario_id, "_replication_", current_replication, ".fasta", sep ="")

# Same for the original fasta file
# Generating saving path for fasta file
list_of_synthetic_population_sequences_original_fasta_current_replication_path   <- paste("R_Output/list_of_synthetic_population_sequences_original_fasta_scenario_id_", scenario_id, "_replication_", current_replication, ".fasta", sep ="")

# Generating saving path for stop codon presence  table
stop_codon_presence_over_time_path   <- paste("R_Output/stop_codon_presence_over_time_scenario_id_", scenario_id, "_replication_", current_replication, ".csv", sep ="")

# Fix 2019.05.04
# Generating saving path for Hamming distance table
hamming_distance_over_time_path   <- paste("R_Output/hamming_distance_over_time_scenario_id_", scenario_id, "_replication_", current_replication, "_short.csv", sep ="")

# Fix 2019.09.10
# Generating saving path for wildtype/mutant status table
wildtype_mutant_status_over_time_path   <- paste("R_Output/wildtype_mutant_status_over_time_scenario_id_", scenario_id, "_replication_", current_replication, "_short.csv", sep ="")




# Fix 2019.06.05
# Generating SHORT versions for saving path for stop codon presence  table
stop_codon_presence_over_time_short_path   <- paste("R_Output/stop_codon_presence_over_time_short_scenario_id_", scenario_id, "_replication_", current_replication, ".csv", sep ="")
# Generating saving path for Hamming distance table
hamming_distance_over_time_short_path   <- paste("R_Output/hamming_distance_over_time_short_scenario_id_", scenario_id, "_replication_", current_replication, "_short.csv", sep ="")
# Fix 2019.09.10.
# Generating saving path for wildtype/mutant status table
wildtype_mutant_status_over_time_short_path   <- paste("R_Output/wildtype_mutant_status_over_time_short_scenario_id_", scenario_id, "_replication_", current_replication, "_short.csv", sep ="")


# Fix 2019.08.26
# To save space on the cluster saving only the string for the sequences after the analysis
# Everything else is NOT saved solely to save the space on the cluster.

# ---- Saving object with the correct namings ----

# Saving object for list of arrays
# eval( parse( text=paste( "save( list_of_synthetic_population_sequences_scenario_id_", scenario_id, "_replication_", current_replication, ", file = list_of_synthetic_population_sequences_current_replication_path )", sep ="") )  )
# Saving object for list of strings
eval( parse( text=paste( "save( list_of_synthetic_population_sequences_strings_scenario_id_", scenario_id, "_replication_", current_replication, ", file = list_of_synthetic_population_sequences_strings_current_replication_path )", sep ="") )  )

# Same for the original array    
# Saving object for list of arrays
# eval( parse( text=paste( "save( list_of_synthetic_population_sequences_original_scenario_id_", scenario_id, "_replication_", current_replication, ", file = list_of_synthetic_population_sequences_original_current_replication_path )", sep ="") )  )
# Saving object for list of strings
# eval( parse( text=paste( "save( list_of_synthetic_population_sequences_original_strings_scenario_id_", scenario_id, "_replication_", current_replication, ", file = list_of_synthetic_population_sequences_original_strings_current_replication_path )", sep ="") )  )



# Saving fasta file
# eval( parse( text=paste( "write.fasta( sequences = list_of_synthetic_population_sequences_strings_scenario_id_", scenario_id, "_replication_", current_replication, ", names = sample_names_synthetic_population, file.out = list_of_synthetic_population_sequences_fasta_current_replication_path )", sep ="") )  )

# Same for the original array    
# Saving fasta file
# eval( parse( text=paste( "write.fasta( sequences = list_of_synthetic_population_sequences_original_strings_scenario_id_", scenario_id, "_replication_", current_replication, ", names = sample_names_synthetic_population, file.out = list_of_synthetic_population_sequences_original_fasta_current_replication_path )", sep ="") )  )



# Fix 2019.06.05.
# Since the files are too big the goal is to save only the small subsets i.e. colums that are incrementin in 100-s.
stop_codon_presence_over_time_short <- stop_codon_presence_over_time[,   list_of_time_slos_to_save ]
hamming_distance_over_time_short   <- hamming_distance_over_time[, list_of_time_slos_to_save ]
# Fix 2019.09.10.
wildtype_mutant_status_over_time_short   <- wildtype_mutant_status_over_time[, list_of_time_slos_to_save ]


# Saving stop codon presence  table
# write.csv( x = stop_codon_presence_over_time, file = stop_codon_presence_over_time_path, row.names = FALSE, quote = FALSE )

# Fix 2019.05.04
# Saving Hamming distance table
# write.csv( x = hamming_distance_over_time,   file = hamming_distance_over_time_path, row.names = FALSE, quote = FALSE )




# Saving stop codon presence  table
write.csv( x = stop_codon_presence_over_time_short,   file = stop_codon_presence_over_time_short_path, row.names = FALSE, quote = FALSE )

# Fix 2019.05.04
# Saving Hamming distance table
write.csv( x = hamming_distance_over_time_short, file = hamming_distance_over_time_short_path,   row.names = FALSE, quote = FALSE )

# Fix 2019.09.10
# Saving wildtype/mutant status table
write.csv( x = wildtype_mutant_status_over_time_short, file = wildtype_mutant_status_over_time_short_path,   row.names = FALSE, quote = FALSE )



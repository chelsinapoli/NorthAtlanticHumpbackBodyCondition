# 2024
# napoli et al. (2024)
# chelsi napoli 

# this code is used to incorporate multiple measurement and altitude errors
# for the morphometric dataset used
# it applies the individual error to all of the individuals
# and assumes that you are working with already wrangled data
# and that you have calculated the altitude error for each drone model & altitude method independently
# using the provided Compare_Accuracy_Script

library(tidyverse)
library(stringr)

# load in your measurements 

# this file should contain your Animal ID, Altitude & Pixel Dimension & Focal Length that your whale was measured using
# the drone model, altitude method (Laser assumes the same laser altimeter model), your width values, your total length
# and then date, DOY, Year, and Foraging Region are useful to carry over throughout the data
meas <- read_csv("Napoli_et_al_Measurements_Pre_Bootstrap.csv")

# # #  Setting Up Your Altitude Error

# For the sake communicating how the code works, I've included how the accuracy code for 2 drone model works below
# You can scale this up to however many drone models x altitude methods you have in your study

# Barometric Altitude 

# run and bring in the results from the Compare Accuracy script
# this script, provided in the Github, allows you to load in your object of known size measurements and compare
# what is expected (the known object's true length) vs observed (what you measured using the drone)
source("Compare_Accuracy_Script_P4P_Baro.R") # this example code is for the Phantom 4 Pro using corrected barometric alts

# The Compare_Accuracy_Script will run, and results in 4 bootstrapped values to find altitude error
# mean altitude error, standard deviation of altitude error, and the 95% confidence interval of altitude error
# pull them into this script, making sure they are numeric 

alt_meanE_P4P <- as.numeric(boot_result_P4P[1]) # mean altitude error for this drone model x alt method
alt_sdE_P4P <- as.numeric(boot_result_P4P[2]) # standard deviation altitude error for this drone model x alt method
alt_025E_P4P <- as.numeric(boot_result_P4P[3]) # 2.5th altitude error for this drone model x alt method
alt_975E_P4P <- as.numeric(boot_result_P4P[4]) # 97.5th altitude error for this drone model x alt method 

boot_means_P4P <- boot_means_P4P # just pulls out the resampled errors for this drone model x alt method so it's not hiding in the R environment 

# Laser Altitude 

# run and bring in the results from the Compare Accuracy script
# this script, provided in the Github, allows you to load in your object of known size measurements and compare
# what is expected (the known object's true length) vs observed (what you measured using the drone)
source("Compare_Accuracy_Script_P4P_Laser.R") # this example code is for the Phantom 4 Pro using laser alts

# The Compare_Accuracy_Script will run, and results in 4 bootstrapped values to find altitude error
# mean altitude error, standard deviation of altitude error, and the 95% confidence interval of altitude error
# pull them into this script, making sure they are numeric 

alt_meanE_P4PL <- as.numeric(boot_result_P4PL[1]) # mean altitude error for this drone model x alt method
alt_sdE_P4PL <- as.numeric(boot_result_P4PL[2]) # standard deviation altitude error for this drone model x alt method
alt_025E_P4PL <- as.numeric(boot_result_P4PL[3]) # 2.5th altitude error for this drone model x alt method
alt_975E_P4PL <- as.numeric(boot_result_P4PL[4]) # 97.5th altitude error for this drone model x alt method 

boot_means_P4PL <- boot_means_P4PL # just pulls out the resampled errors for this drone model x alt method so it's not hiding in the R environment 

# Convert All Measurements to Pixels 
# This allows you scale your measurements with an error-adjusted altitude

# this converts your total length (col 5) and all width measurements (6 to 24) to pixels
allpixels_og <- meas %>% mutate(across(5:24, ~ . * `Focal Length` / (Altitude * `Pixel Dimension`)))

# # # Setting Up Inter-Image Error

# Here is where you're going to load in the data that you have that includes the additional measurements
# of the same individuals 
# This is to better account for and propagate error from image quality

# Load in additional measurements

addi_join <- read.csv("Napoli_et_al_Additional_Measurements.csv") # read in the CSV file 

addi_join <- addi_join %>% dplyr::select(-X) # get rid of the additional column added in
meas <- meas %>% dplyr::select(-`...1`) # get rid of the additional column added in 

colnames(addi_join) = colnames(meas) # make sure that the column names are the same and represent the same columns

# this adds an identifier for the source of the data, if it was the best quality image (OG = original) or additional (Addi)
addi_join <- addi_join %>% mutate(SRC = "Addi") %>% mutate(srcID = paste0(Animal_ID, "_", SRC))
alldata_og <- meas %>% mutate(SRC = "OG") %>% mutate(srcID = paste0(Animal_ID, "_", SRC))

# clean up any binding conflicts between different types of data; the date isn't needed in Date format for this so it's easiest to change it to a character
alldata_og$Date <- as.character(alldata_og$Date)

alldata_comp <- bind_rows(alldata_og, addi_join) %>% arrange(Animal_ID) # this binds the rows and arranges by Animal ID 

alldata_dupesonly <- alldata_comp %>% group_by(Animal_ID) %>% filter( n() > 1 ) # this takes out any animals that do not have a duplicate measurement

# # ! ! ! ! ! THIS IS AN INCREDIBLY IMPORTANT NOTE. THE BELOW CODE ASSIGNS THE ORIGINAL, BEST MEASUREMENT AS THE "MEAN" 
# # IT IS NOT! TAKING! THE! MEAN! VALUE! 
# # THIS CODE IS SIMPLY COPYING THE OG VALUE AND ASSIGNING IT AS A "MEAN" VALUE BUT IT'S TAKING THE "MEAN" OF ONLY 1 VALUE, THE ORIGINAL MEASUREMENT.
# # THIS IS WHERE THE CODE IS NOT LABELED INTUITIVELY (AS "MEAN" INSTEAD OF "ORIGINAL" BUT THIS CODE WAS ADAPTED AND
# # ADJUSTING ALL THE NAMES AND FUNCTIONS BROKE THE CODE # I'm sorry! I know this isn't intuitive. 
# # IMPORTANTLY, THE MEAN VALUE HERE IS NOT THE MEAN OF DUPLICATE MEASUERMENT, BUT THE ORIGINAL BEST MEASUREMENT, AS DESCRIBED ACCURATELY IN THE METHODS OF NAPOLI ET AL. 2024
# # Apologies for any confusion on this, but it's also very easy to verify this by using the "View(alldata_meanTL)" function and seeing that the
# # "mean" values are indeed the original measurement instead 

# # Just to iterate from the methods of Napoli et al (2024): The differences in total length and width measurements between the best scoring image and additional
# # images of each whale were calculated (Figure 2). Measurement difference was sampled with replacement 1000 times to form an error distribution for the total length 
# # and each 5% width for each 747 drone model, except the DJI Mavic 3, DJI Inspire 1, and DJI Inspire 2 in which multiple measurements were not available or sufficient to form an error distribution (n = 0, 1, and 2, respectively).

alldata_meanTL <-  alldata_dupesonly %>% 
  group_by(Animal_ID) %>% # grouping by animal ID
  mutate(across(TL:`95% Width`, list(mean = ~mean(.[SRC == "OG"], na.rm = TRUE)))) # once again, not taking the mean of duplicates but copying the OG measurement

# Set up the Bootstrap Function 

bootstrap_func_meas <- function(mean_meas, measured_meas, iteration_number){
  error <- mean_meas - measured_meas # error between OG measurement and measured measurement 
  bootstrap_means <- numeric(iteration_number) # number of iterations
  
  for(i in 1:iteration_number){
    # Sample with Replacement
    bootstrap_sample <- sample(error, replace = TRUE)
    # Calculate Mean Error For BootStrap Sample
    bootstrap_means[i] <- mean(bootstrap_sample)
    
  }
  
  # Bootstrapped Mean Error
  bootstrapped_mean_error <- mean(bootstrap_means) 
  # Calculate Confidence Interval
  bootstrapped_CI <- quantile(bootstrap_means, c(0.025, 0.975))
  bootstrapped_sd <- sd(bootstrap_means)
  
  return(c(bootstrapped_mean_error, bootstrapped_CI, bootstrapped_sd)) # return the bootstrapped mean error, bootstrapped confidence interval, and bootstrapped standard deviation
  
}

# # # 

# Set up A Loop for Each Drone Model At Each Width 

# SAME THING HERE. "MEANS" REFERS TO THE ORIGINAL MEASUREMENTS, NOT AN ACTUAL MEAN. IT IS TECHNICALLY THE MEAN OF THE ORIGINAL MEASUREMENT, WHICH IS ITSELF.
pairsnames_meas <- colnames(alldata_meanTL)[4:23] # pulling out the measurements for all columns from Total_Length to 95% Width # double check your col numbers
pairsnames_means <- colnames(alldata_meanTL)[34:53] # pulling out the OG measurements for all columns from Total_Length to 95% Width # double check your col numbers

interimage_boot_list <- list() # set up a list for the loop to land into

unique_drone_models <- unique(alldata_meanTL$Drone_Model) # what are the unique drone models # aka what are the unique camera systems being used

# This Loop Generates the Inter-Image Resampling Results 
for (j in 1:length(unique_drone_models)){ # for each drone model 
  
  dronefilt_meanTL <- alldata_meanTL %>% filter(Drone_Model == unique_drone_models[j])  # filter data for just one drone model in the vector
  
  interimage_boot_results <- data.frame() # create a dataframe for the results of that drone model 
  for (i in 1:length(pairsnames_meas)){ # for each width 
    
    mean_meas <- dronefilt_meanTL %>% group_by(Animal_ID) %>% dplyr::pull(pairsnames_means[i]) # pulls the OG measurements 
    measured_meas <- dronefilt_meanTL %>% group_by(Animal_ID) %>% dplyr::pull(pairsnames_meas[i]) # pull the measured measurements
    
    boot_result <- bootstrap_func_meas(mean_meas = mean_meas, measured_meas = measured_meas, iteration_number = iteration_number) # bootstrap function
    
    error_mean <- as.numeric(boot_result[1]) # mean resampled error as numeric
    error025 <- as.numeric(boot_result[2]) # resampled error 2.5%
    error975 <- as.numeric(boot_result[3]) # resampled error 97.5% 
    error_sd <- as.numeric(boot_result[4]) # resampled standard deviation
    
    row <- c(pairsnames_meas[i], error_mean, error_sd, error025, error975) # a data row for each width, the mean, sd, and CI
    interimage_boot_results <- rbind(row, interimage_boot_results) # put into the results 
  }
  
  colnames(interimage_boot_results) <- c("Width", "Mean", "SD", "Error_025", "Error_975") # reformat the data
  
  convert_to_numeric <- function(df, columns) {
    df <- df %>%
      mutate_at(vars(columns), as.numeric)
    return(df)
  } # convert to numeric just to be extra sure it's in numeric format 
  
  interimage_boot_results <- convert_to_numeric(interimage_boot_results, columns = c("Mean", "SD", "Error_025", "Error_975")) # convert 
  
  interimage_boot_results[2:5] <- round(interimage_boot_results[2:5], 9) # round 
   
  interimage_boot_list[[j]] <- interimage_boot_results # save the results for every drone model into a list 
  
}


# # For the Mavic 3, Inspire 1, and Inspire 2 there isn't enough repeat measurement data to actually form a distribution that is valid, 
# so we used the measurement error that is pooled from the other images
# note that the above code does actually produce distribution information for Inspire 1 and Inspire 1 but that error distribution is formed from
# minimal measurements, so it's technically computed but would not be scientifically sound to use due to such a small sample size 

# So here's that same code, with the same functionality, just pooling the data as opposed to separating it out into unique drone models
# hence there is no nested loop at the j-level 

interimage_boot_results_pooled <- data.frame()
for (i in 1:length(pairsnames_meas)){
  
  mean_meas <- alldata_meanTL %>% group_by(Animal_ID) %>% dplyr::pull(pairsnames_means[i])
  measured_meas <- alldata_meanTL %>% group_by(Animal_ID) %>% dplyr::pull(pairsnames_meas[i])
  
  boot_result <- bootstrap_func_meas(mean_meas = mean_meas, measured_meas = measured_meas, iteration_number = iteration_number)
  
  error_mean <- as.numeric(boot_result[1])
  error025 <- as.numeric(boot_result[2])
  error975 <- as.numeric(boot_result[3])
  error_sd <- as.numeric(boot_result[4])
  
  row <- c(pairsnames_meas[i], error_mean, error_sd, error025, error975)
  interimage_boot_results_pooled <- rbind(row, interimage_boot_results_pooled)
}

colnames(interimage_boot_results_pooled) <- c("Width", "Mean", "SD", "Error_025", "Error_975")

convert_to_numeric <- function(df, columns) {
  df <- df %>%
    mutate_at(vars(columns), as.numeric)
  return(df)
}

interimage_boot_results_pooled <- convert_to_numeric(interimage_boot_results_pooled, columns = c("Mean", "SD", "Error_025", "Error_975"))

interimage_boot_results_pooled[2:5] <- round(interimage_boot_results_pooled[2:5], 9)

# Up until this bit, it's the same as above 

interimage_boot_list[[7]] <- interimage_boot_results_pooled # here is where we add the pooled results to the unique drone model results list above

names(interimage_boot_list) <- c(unique_drone_models, "Mavic 3") # Here, we name the interimage boot list above using the corresponding drone models, and assign the pooled results (7) as the Mavic

interimage_boot_list[["Inspire 1"]] <- interimage_boot_results_pooled # Here, we replace the Inspire 1 that was generated using a too small sample size with the pooled data
interimage_boot_list[["Inspire 2"]] <- interimage_boot_results_pooled # Here, we replace the Inspire 2 that was generated using a too small sample size with the pooled data 

# # # # 

# Okay, so here's the big error distribution generating loop
# In summary, this loop generates the 1000 populations worth of humpback whale measurements that incorporate altitude error from each drone model x alt method combo
# and it incorporates the inter-image error for each width for the corresponding drone model used
# this loop it going to take your computer a couple hours to run # it's a long time
# after you create this list, you will run all of the 1000 x n measurements (in this case, 256 x 1000 = 256,000) whales through the 3D model
# that step, with 256k whales, usually takes about 24 hours to run using the blender code from Hirtle et al. 
# note that these 256k whales are THEN randomly sampled into 1000 populations again (next code set available on the Github)

Bootstrap_Error_List <- list() # create a list for storage 

for(m in 1:1000){ # this is the big loop that will do this 1000 times 
  
  Error_Adj_Meas <- data.frame() # data frame to store the error adjusted measurements
  for (i in 1:length(allpixels_og$Animal_ID)){ # using the original measurements in pixels! important to use pixels here because you need to re-scale your pixels using the new altitudes 
    
    animal_ID <- allpixels_og$Animal_ID # animal ID 
    ref_row <- allpixels_og[i, ] # the reference row
    drone_error <- ref_row %>%
      mutate(Drone_Error = case_when(
        Drone_Model == "Phantom 4 Std" & Alt_Method == "Barometric" ~ sample(boot_means_GOMP4S, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        Drone_Model == "Phantom 4 Pro" & Alt_Method == "Barometric" ~ sample(boot_means_P4P, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        Drone_Model == "Matrice 200" & Alt_Method == "Barometric" ~ sample(boot_means_GL_MAT200, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        Drone_Model == "Inspire 1" & Alt_Method == "Barometric" ~ sample(boot_means_IN1, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        Drone_Model == "Phantom 4 Pro" & Alt_Method == "Laser" ~ sample(boot_means_P4PL, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        Drone_Model == "Mavic 3" & Alt_Method == "Laser" ~ sample(boot_means_MAVIC, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        Drone_Model == "Matrice 210" & Alt_Method == "Laser" ~ sample(boot_means_MAT210, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        Drone_Model == "Inspire 2" & Alt_Method == "Laser" ~ sample(boot_means_MAT210, 1), # if this is the corresponding drone x alt combination to that reference row, this samples 1 value from the bootstrapped differences between true value and measured value of an object of known size 
        TRUE ~ as.numeric(0))) %>% # this is here only as a fail-safe for the code, that if there is some other combination, there's no error.
      dplyr::pull(Drone_Error) # pull out the corresponding sampled altitude error from each drone model x altimeter method combo 
    
    alt_with_error <- ref_row$Altitude+drone_error # this adds the altitude error from the model x alt method combo to the original altitude 
    
    ref_row$alt_with_error <- alt_with_error # adds it to the reference row 
    
    ALT_ref_row <- ref_row %>%
      mutate(across(6:25, ~ . * ((alt_with_error * `Pixel Dimension`)/`Focal Length`))) # this adjusts all of the measurements to the new altitude that has altitude error incorporated into it
    
    # This is the end of incorporating altitude error
    
    # Now to incorporate specific measurement error of the drone
    
    drone <- ALT_ref_row$Drone_Model # this is the drone model used for the reference row that we're working with
    interimage_boot_df <- interimage_boot_list[[drone]] # this pulls out the mean error for each width for the drone model 
    
    interimage_error <- data.frame() # data frame to store 
    for(j in 1:length(interimage_boot_df$Width)){ # this loop looks at each width 
      width_row <- interimage_boot_df[j,] # this is the row with the widths in it
      width <- width_row$Width # the width of choice 
      error <- rnorm(1, mean = width_row$Mean, sd = width_row$SD) # sampling the width error from the distribution with the drone's parameters
      
      row <- c(width, error) # forming a new row
      interimage_error <- rbind(row, interimage_error) # binding it 
    }
    
    colnames(interimage_error) <- c("Width", "Interimage_Error") # creating a labeled data frame for the width and corresponding error 
    
    addi_error <- interimage_error %>%
      mutate(Interimage_Error = as.numeric(Interimage_Error)) %>% # make sure it's numeric 
      pivot_wider(names_from = Width, values_from = Interimage_Error) # pivot from the width to the interimage error 
    
    adj_error_df <- data.frame() # make new data frame to store it 
    
    # for each width in the reference row that has been scaled with the error-adjusted-new-altitude, take the width and add the error to it
    for(k in 1:length(colnames(addi_error))){ 
      widths <- colnames(addi_error)
      width <- widths[k]
      adj_error <- ((ALT_ref_row %>% dplyr::pull(width))) + addi_error %>% dplyr::pull(width)
      temp <- data.frame(Name = width, Adj_Error = adj_error)
      adj_error_df <- rbind(temp, adj_error_df)
    }
    
    adj_error_join <- adj_error_df %>% pivot_wider(names_from = Name, values_from = Adj_Error) %>% 
      dplyr::select(colnames(ref_row)[6:25]) # reformat and replace with error adjusted row
    
    adj_ref_row <- sjmisc::replace_columns(ref_row, adj_error_join) # reformat and replace with error adjusted row 
    
    Error_Adj_Meas <- rbind(adj_ref_row, Error_Adj_Meas) # save each error-adjusted-row to the data frame 
  }
  
  Error_Adj_Meas$Pop_ID <- paste0("PopID","_", m) # save each iteration with a unique PopID identfifier # NOTE: this PopID is NOT the same as the population that will get analyzed # This PopID is to ensure that the resampled populations are randomly sampled and not just iterative # this ensures 2 levels of random sampling
  
  Bootstrap_Error_List[[m]] <- Error_Adj_Meas # save in the list
  
}

Bootstrap_MegaDF <- bind_rows(Bootstrap_Error_List) # bind all of the 256k whales that have had resampled error incorporated into a huge data frame that can be run through the Hirtle et al. Blender model 


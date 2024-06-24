# 2024
# napoli et al. (2024)
# chelsi napoli 

# this code is used to compare the known length of an object with the measured length of an object using a drone
# this code can be adapted for every drone model x altimeter method combination in your study
# for napoli et al (2024) there were 8 versions of this code that only differed in the folders they linked to and the true length of the measurement
# this code assumes that you used morphometrix (https://github.com/wingtorres/morphometrix) to make your measurements of a known object
# Torres, W.I., and Bierlich, K.C (2020). MorphoMetriX: a photogrammetric measurement GUI for morphometric analysis of megafauna.. Journal of Open Source Software, 4(44), 1825. https://doi.org/10.21105/joss.01825

library(tidyverse)
library(stringr)

set.seed(03282020)

# Collect All of the CSV Measurement files that correspond to the drone model and altimeter method used
# This assumes that your measurement files of the known object are separated into folders 

meas_csvs <- list.files(path = "./P4P_Laser/Measurements", pattern = ".csv", recursive = TRUE, full.names = TRUE)

# Function to Munge Morphs 
munge_morphs <- function(x){
  d_meta <- read_csv(file = x, n_max = 9, col_names = FALSE)
  d_meta <- spread(d_meta, X1, X2) %>% mutate(`Whale ID` = `Image ID`) %>% dplyr::select(`Image Path`, `Whale ID`, Altitude, `Focal Length`, `Pixel Dimension`, TL) %>% dplyr::rename("Pixel_Dimension" = "Pixel Dimension", "Focal_Length" = "Focal Length")
  
  return(d_meta)
}

# This function creates a data frame that binds all of the Morphometrix outputs together, and forms a dataframe with column for: 
# Image Path (where your known object image lives), Whale ID (However you labeled your photo in Morphometrix at measurement time), Altitude, Focal Length, Pixel Dimension, and Total Length (TL, or however you called your length in Morphometrix)
meas_meas <- bind_rows(lapply(X = meas_csvs, FUN = munge_morphs))

# Check your data before using! # Sometimes the munge function creates extra lines that are all NA, this function gets rid of them
# Double check to make sure that you're not removing images that should have corresponding measurements that are missing! that is NOT the purpose of this line of code 
meas_meas <- meas_meas %>% filter(!is.na(TL)) 

meas_meas$BaseName <- basename(meas_meas$`Image Path`) # extracts just the filename, gets rid of the rest of the file path 

compare <- meas_meas # creates a copy of the data frame, not necessary but is nice for any potential debugging 

TL_Correct <- 1 # this is the true length of your known object in the same unit as Morphometrix (m)

compare$TL_Correct <- as.numeric(TL_Correct) # make sure it's numeric 
compare$TL <- as.numeric(compare$TL) # make sure it's numeric 
compare <- compare %>% filter(complete.cases(.)) # once again, a good check in case of any lines that are NAs from the munge function

# Differences 

compare_diff <- compare %>% mutate(Meas_Diff = TL - TL_Correct) # this computes your measurement different 

compare_diff$Focal_Length <- as.numeric(compare_diff$Focal_Length) # make sure your values are numeric 
compare_diff$Pixel_Dimension <- as.numeric(compare_diff$Pixel_Dimension) # make sure your values are numeric 
compare_diff$Altitude <- as.numeric(compare_diff$Altitude) # make sure your values are numeric 

compare_diff <- compare_diff %>% 
  mutate(Pixels = (Focal_Length*TL)/(Altitude*Pixel_Dimension)) # this calculates the length in number of pixels # essentially back calculates Morphometrix's output into pixels

compare_corrected <- compare_diff %>%
  mutate(Correct_Alt = (TL_Correct * Focal_Length)/(Pixels*Pixel_Dimension)) # finds the true altitude of the object needed to measure the pixels and result in the object's true known length

compare_alt_diff <- compare_corrected %>%
  mutate(Alt_Diff_Corrected = Correct_Alt - Altitude) # finds the difference in the corrected altitude and the altitude from the drone model x laser altimeter combination 

# Now we create two functions: 

# First function returns a summary of the bootstrapping: returns the mean error, mean confidence interval, and mean standard deviation
# This is important because this is what we will use to report the bootstrapping errors in a summarized way
bootstrap_func <- function(og, corrected, iteration_number){
  error <- og - corrected # difference between original and corrected value 
  bootstrap_means <- numeric(iteration_number) # number of times to be resampled 
  
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
  bootstrapped_SD <- sd(bootstrap_means)
  
  return(c(bootstrapped_mean_error, bootstrapped_SD, bootstrapped_CI)) # This is the key part of this function, and what separates it from the first bootstrap function
  
}

# Second function returns the actual differences between the original and corrected 1000 times - this is the raw, un-summarized resampling of the errors
# This is important because this is what we will sample from directly when we incorporate and propagnate the error into the measurements
bootstrap_means_func <- function(og, corrected, iteration_number){
  error <- og - corrected
  bootstrap_means <- numeric(iteration_number)
  
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
  bootstrapped_SD <- sd(bootstrap_means)
  
  return(c(bootstrap_means)) # This is the key part of this function, and what separates it from the first bootstrap function
  
  
}

iteration_number <- 1000 # set your iteration number 

# This is the summarized stats of the bootstrapped error for this drone model x altimeter method 
boot_result_P4PL <- bootstrap_func(og = compare_alt_diff$Altitude, corrected = compare_alt_diff$Correct_Alt, iteration_number = iteration_number)

# This is the resampled errors (boot_means_P4PL should have a length(boot_means_P4PL) the same as your iteration number (1000))
# This is what we will pull from in the error-incorporating-for-individuals script (Napoli_et_al_Bootstrapping_Individuals.R)
boot_means_P4PL <- bootstrap_means_func(og = compare_alt_diff$Altitude, corrected = compare_alt_diff$Correct_Alt, iteration_number = iteration_number)


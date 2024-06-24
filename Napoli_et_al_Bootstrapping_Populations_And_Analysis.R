# # 2024
# napoli et al. (2024)
# chelsi napoli 

# this code is to bootstrap populations that have incorporated altitude and interimage error into measurements
# compute a Body Condition Index for each whale in each population, and analyze the data using ANCOVAs and SMAs

# this code is a bit of a monster ! sorry for it's length! 
# also, it's pretty hyper-specific for napoli et al. 2024, especially the loops and functions that pull out data according to
# region which, as a reminder, is sorted alphabetically as a default in R
# so I've tried to avoid hard-coding things as much as I can so that other people can use the code "straight out of the box"
# generally there are notes where things are hardcoded because of the specifics of the project as opposed to what should actually be hardcoded


# Setup -------------------------------------------------------------------

set.seed(03282020)

# to be honest, I don't think you need all these packages for this code but as I've adapted the code I can't remember which packages I need and which ones I don't
# so they've all stayed thoughout all of the versions
library(tidyverse)
library(lubridate)
library(car)
library(smatr)
library(dplyr)
library(npreg)
library(stringr)
library(ggformula)
library(multcomp)
library(readxl)
library(ggsignif)
library(plyr)
library(coin)
library(FSA)
library(bootstrap)
library(sm)
library(lsr)
library(lme4)
library(readr)

# The code assumes that you've generated your error-incorporated measurements using the Napoli_et_al_Bootstrapping_Individuals code
# and that you've run all of your models through Hirtle et al. Blender code to generate 3D models and segment & whole body volume estimates

# It assumes that you've added up the volume segments of the whale at each 5% slice to compute a whole body volume 

ind_vols_total <- read.csv("Your_3D_Model_Output.csv")

# If you haven't added up the volume segments of the whale at each 5% slice, here's some code to do that
ind_vols <- ind_vols_total %>%
            dplyr::mutate(Body_Volume = rowSums(across(`0`:`34`))) %>% # add up the rows for total Body Volume 
            dplyr::select(-c(`0`:`34`)) %>% # remove the columns with individual segment volumes 
            dplyr::rename("Total_Length" = "TL") # rename Total Length 


# Assemble Populations ----------------------------------------------------

# Assemble Populations 

# This code forms 1000 populations of whales by randomly selecting error-incorporated measurements that correspond to each unique animal ID
# So this creates 1000 population of whales all with 256 whales. Each of the 256 whales are reflective of the best quality images that
# have both interimage and altitude error incorporated into the measurements, and body volume estimates that were derived from each drone and altimeter model combination

pop_boot_list <- list() # Empty list to store the results 
pop_boot_df <- data.frame() # Empty data frame to store the results 
for(i in 1:1000){
  all_Animal_ID <- unique(ind_vols$Animal_ID) # all of the unique Animal ID 
  pop_boot_temp <- data.frame() # temp data frame 
  for (j in 1:length(all_Animal_ID)){
    extract <- ind_vols %>% filter(Animal_ID == all_Animal_ID[j]) # Filter all of the rows with the unique Animal ID 
    row_sample <- sample_n(extract, size = 1, replace = TRUE) # Sample randomly 1 of the unique Animal ID rows 
    pop_boot_temp <- rbind(row_sample, pop_boot_temp) # Bind that into the data frame 
  }
  pop_boot_df <- pop_boot_temp # each population of 256 stored here
  pop_boot_df$Population_Number <- i # each population number assigned 
  pop_boot_list[[i]] <- pop_boot_df # each of the 1000 populations with pop number stored into a list
}


# This list takes a long time to generate (1 hr)
# I've loaded in the .RDS below of the list available on Github that was generated using the same set.seed(03282020) for Napoli et al 2024
# A .csv version of the same list with rows bound together is available on Github, though this code utitizes lapply for many of the analyses so
# having it not in list format is relatively useless, and I suggest converting the .csv back to list by Population_Number if necessary 

pop_boot_list <- readRDS("Napoli_et_al_Bootstrapped_Populations_RList.RDS")








# Set Up Functions --------------------------------------------------------

# The analysis for all of the bootstrapped data really hinges on using the lapply function
# This is for two reasons: 1) All of the 1000 populations are in a list already 2) The analysis is identical for each population

# The first part of this code is really just setting up the functions so that they can later be used
# using the lapply function 

# Convert to Numeric - Takes a dataframe and the columns you want to convert to numeric and does that
convert_to_numeric <- function(df, columns) {
  df <- df %>%
    mutate_at(vars(columns), as.numeric)
  return(df)
}

# Calculate Body Condition Index
# Assigns age class based on Bierlich's 11.47 m unknown sex value
# Uses BCI equation
# Returns dataframe with BCI calculated for each individual in the population 
BCI_func <- function(df){
  
  temp <- df
  temp <- temp %>% mutate(age_class = ifelse(Total_Length > 11.47, "Adult", "Juvenile")) # cut off of adult being 11.47 meters as used in Bierlich. (2022)
  temp$Year <- as.factor(as.character(temp$Year))
  
  bci_equation <- lm(log(temp$Body_Volume) ~ log(temp$Total_Length)) # equation to find expected body volume from Christiansen et al (2020)
  
  temp$logfit <- fitted(bci_equation) # logfit 
  temp$fit <- exp(temp$logfit) # expected value 
  temp$BCI <- (temp$Body_Volume-temp$fit)/temp$fit # equation for BCI 
  
  return(temp)
}

# Run Age Class T-Tests
# Extracts p-values, t-values, degrees of freedrom
# Creates summary data frame with t-test stats
ageclass_ttest_func <- function(df){
  temp <- df
  nyb <- temp %>% filter(Foraging_Region == "New York Bight")
  gom <- temp %>% filter(Foraging_Region == "Gulf of Maine")
  ice <- temp %>% filter(Foraging_Region == "Iceland")
  gl <- temp %>% filter(Foraging_Region == "Greenland")
  ttest_nyb <- t.test(BCI ~ rep_class, data = nyb)
  ttest_gom <- t.test(BCI ~ rep_class, data = gom)
  #ttest_ice <- t.test(BCI ~ rep_class, data = ice) # removed because not enough sample size
  ttest_gl <- t.test(BCI ~ rep_class, data = gl)
  
  nyb_pvalue <- ttest_nyb$p.value
  nyb_tvalue <- as.numeric(ttest_nyb$statistic)
  nyb_df <- as.numeric(ttest_nyb$parameter)
  gom_pvalue <- ttest_gom$p.value
  gom_tvalue <- as.numeric(ttest_gom$statistic)
  gom_df <- as.numeric(ttest_gom$parameter)
  # ice_pvalue <- ttest_ice$p.value
  # ice_tvalue <- as.numeric(ttest_ice$statistic)
  # ice_df <- as.numeric(ttest_ice$parameter)
  gl_pvalue <- ttest_gl$p.value
  gl_tvalue <- as.numeric(ttest_gl$statistic)
  gl_df <- as.numeric(ttest_gl$parameter)
  
  ttest_results_t <- as.data.frame(rbind(nyb_pvalue, nyb_tvalue, nyb_df,
                                         gom_pvalue, gom_tvalue, gom_df,
                                         #ice_pvalue, ice_tvalue, ice_df,
                                         gl_pvalue, gl_tvalue, gl_df))
  ttest_results <- data.table::transpose(ttest_results_t)
  colnames(ttest_results) <- rownames(ttest_results_t)
  
  return(ttest_results)
  
}

# Run Age Class ANCOVA
# Models BCI as a function of Age Class with Foraging Region as a random effect 
ancova_agemodel_func <- function(df){
  temp <- df
  temp$Foraging_Region <- as.factor(temp$Foraging_Region)
  temp$DOY <- as.numeric(temp$DOY)
  agemodel <- aov(BCI ~ age_class + Error(Foraging_Region), data = temp)
  return(agemodel)
}

# This function extracts the statistics from the Age Class ANCOVA Results
# Makes a summary table and returns them
extract_ancova_agemodel_func <- function(list){
  ancova_agemodel <- list
  ancova_agemodel_summary <- summary(ancova_agemodel)
  error_within <- unlist(ancova_agemodel_summary$`Error: Within`)
  error_foragingregion <- unlist(ancova_agemodel_summary$`Error: Foraging_Region`)
  error_within_DF1 <- as.numeric(error_within[1])
  error_within_DF2 <- as.numeric(error_within[2])
  error_within_SumSq1 <- as.numeric(error_within[3])
  error_within_SumSq2 <- as.numeric(error_within[4])
  error_within_MeanSq1 <- as.numeric(error_within[5])
  error_within_MeanSq2 <- as.numeric(error_within[6])
  error_within_FValue1 <- as.numeric(error_within[7])
  error_within_FValue2 <- as.numeric(error_within[8])
  error_within_Pvalue1 <- as.numeric(error_within[9])
  error_within_Pvalue2 <- as.numeric(error_within[10])
  error_foragingregion_DF1 <- as.numeric(error_foragingregion[1])
  error_foragingregion_DF2 <- as.numeric(error_foragingregion[2])
  error_foragingregion_SumSq1 <- as.numeric(error_foragingregion[3])
  error_foragingregion_SumSq2 <- as.numeric(error_foragingregion[4])
  error_foragingregion_MeanSq1 <- as.numeric(error_foragingregion[5])
  error_foragingregion_MeanSq2 <- as.numeric(error_foragingregion[6])
  error_foragingregion_FValue1 <- as.numeric(error_foragingregion[7])
  error_foragingregion_FValue2 <- as.numeric(error_foragingregion[8])
  error_foragingregion_Pvalue1 <- as.numeric(error_foragingregion[9])
  error_foragingregion_Pvalue2 <- as.numeric(error_foragingregion[10])
  
  ancova_results_t <- as.data.frame(rbind(error_within_DF1,error_within_DF2,
                                          error_within_SumSq1,error_within_SumSq2,
                                          error_within_MeanSq1,error_within_MeanSq2,
                                          error_within_FValue1,error_within_FValue2,
                                          error_within_Pvalue1,error_within_Pvalue2,
                                          error_foragingregion_DF1,error_foragingregion_DF2,
                                          error_foragingregion_SumSq1,error_foragingregion_SumSq2,
                                          error_foragingregion_MeanSq1,error_foragingregion_MeanSq2,
                                          error_foragingregion_FValue1,error_foragingregion_FValue2,
                                          error_foragingregion_Pvalue1,error_foragingregion_Pvalue2))
  ancova_results <- data.table::transpose(ancova_results_t)
  colnames(ancova_results) <- rownames(ancova_results_t)
  
  return(ancova_results)
  
}

# This function runs the main ANCOVA that looks at the relationship between BCI and Foraging Region with covariates for DOY and Year
# Also runs the Tukey Post Hoc Test
# And also runs the effect size 
# Returns all three in a list
ancova_func <- function(df){
  temp <- df
  temp$Foraging_Region <- as.factor(temp$Foraging_Region)
  temp$DOY <- as.numeric(temp$DOY)
  ancova_resid <- aov(BCI ~ Foraging_Region + DOY + Year, data=temp)
  ancova_results <- Anova(ancova_resid, type = "III")
  
  #tukey post hoc 
  postHocs_resid <- glht(ancova_resid, linfct = mcp(Foraging_Region = "Tukey"))
  postHocs_summary <- summary(postHocs_resid)
  
  #effect size 
  etaSq <- etaSquared(ancova_resid)
  
  return(list(ancova_results, postHocs_summary, etaSq))
  
}

# This function extracts the statistics for the ANCOVA between BCI ~ Foraging Region + DOY + Year
# Creates a table and summarizes the statistics
# Returns the statistics 
extract_ancova_func <- function(list){
  ancova_table <- list[[1]]
  intercept_sum_sq <- ancova_table$`Sum Sq`[1]
  intercept_df <- ancova_table$Df[1]
  intercept_fvalue <- ancova_table$`F value`[1]
  intercept_pval <- ancova_table$`Pr(>F)`[1]
  Foraging_Region_sum_sq <- ancova_table$`Sum Sq`[2]
  Foraging_Region_df <- ancova_table$Df[2]
  Foraging_Region_fvalue <- ancova_table$`F value`[2]
  Foraging_Region_pval <- ancova_table$`Pr(>F)`[2]
  DOY_sum_sq <- ancova_table$`Sum Sq`[3]
  DOY_df <- ancova_table$Df[3]
  DOY_fvalue <- ancova_table$`F value`[3]
  DOY_pval <- ancova_table$`Pr(>F)`[3]
  Year_sum_sq <- ancova_table$`Sum Sq`[4]
  Year_df <- ancova_table$Df[4]
  Year_fvalue <- ancova_table$`F value`[4]
  Year_pval <- ancova_table$`Pr(>F)`[4]
  
  ancova_results_t <- as.data.frame(rbind(intercept_sum_sq, intercept_df, intercept_fvalue, intercept_pval,
                                          Foraging_Region_sum_sq, Foraging_Region_df, Foraging_Region_fvalue, Foraging_Region_pval,
                                          DOY_sum_sq, DOY_df, DOY_fvalue, DOY_pval,
                                          Year_sum_sq, Year_df, Year_fvalue, Year_pval))
  ancova_results <- data.table::transpose(ancova_results_t)
  colnames(ancova_results) <- rownames(ancova_results_t)
  
  return(ancova_results)
  
}

# This function extracts the statistics from the Tukey test from BCI ~ Foraging Region + DOY + Year
# Creates a table and summarizes the statistics
# Returns the statistics 
# Sorry that this is so hard coded based on foraging region
# But essentially it's in alphabetical order - with GOM [1], GReenland [2], Iceland [3], New York Bight [4]
extract_tukey_func <- function(list){
  tukey_table <- list[[2]]
  gom_greenland_estimate <- as.numeric(tukey_table$test$coefficients[1])
  gom_greenland_stderror <- as.numeric(tukey_table$test$sigma[1])
  gom_greenland_tvalue <- as.numeric(tukey_table$test$tstat[1])
  gom_greenland_pvalue <- as.numeric(tukey_table$test$pvalues[1])
  iceland_greenland_estimate <- as.numeric(tukey_table$test$coefficients[2])
  iceland_greenland_stderror <- as.numeric(tukey_table$test$sigma[2])
  iceland_greenland_tvalue <- as.numeric(tukey_table$test$tstat[2])
  iceland_greenland_pvalue <- as.numeric(tukey_table$test$pvalues[2])
  nyb_greenland_estimate <- as.numeric(tukey_table$test$coefficients[3])
  nyb_greenland_stderror <-as.numeric(tukey_table$test$sigma[3])
  nyb_greenland_tvalue <- as.numeric(tukey_table$test$tstat[3])
  nyb_greenland_pvalue <- as.numeric(tukey_table$test$pvalues[3])
  iceland_gom_estimate <- as.numeric(tukey_table$test$coefficients[4])
  iceland_gom_stderror <- as.numeric(tukey_table$test$sigma[4])
  iceland_gom_tvalue <- as.numeric(tukey_table$test$tstat[4])
  iceland_gom_pvalue <- as.numeric(tukey_table$test$pvalues[4])
  nyb_gom_estimate <- as.numeric(tukey_table$test$coefficients[5])
  nyb_gom_stderror <- as.numeric(tukey_table$test$sigma[5])
  nyb_gom_tvalue <- as.numeric(tukey_table$test$tstat[5])
  nyb_gom_pvalue <- as.numeric(tukey_table$test$pvalues[5])
  nyb_iceland_estimate <- as.numeric(tukey_table$test$coefficients[6])
  nyb_iceland_stderror <- as.numeric(tukey_table$test$sigma[6])
  nyb_iceland_tvalue <- as.numeric(tukey_table$test$tstat[6])
  nyb_iceland_pvalue <- as.numeric(tukey_table$test$pvalues[6])
  tukey_results_t <- as.data.frame(rbind(     gom_greenland_estimate, gom_greenland_stderror, gom_greenland_tvalue, gom_greenland_pvalue,
                                              iceland_greenland_estimate, iceland_greenland_stderror, iceland_greenland_tvalue, iceland_greenland_pvalue,
                                              nyb_greenland_estimate, nyb_greenland_stderror, nyb_greenland_tvalue,nyb_greenland_pvalue,
                                              iceland_gom_estimate,iceland_gom_stderror, iceland_gom_tvalue,iceland_gom_pvalue,
                                              nyb_gom_estimate,nyb_gom_stderror, nyb_gom_tvalue,nyb_gom_pvalue,
                                              nyb_iceland_estimate,nyb_iceland_stderror,nyb_iceland_tvalue,nyb_iceland_pvalue))
  tukey_results <- data.table::transpose(tukey_results_t)
  colnames(tukey_results) <- rownames(tukey_results_t)
  
  return(tukey_results)
  
}

# This function extracts the statistics from the EtaSquare from BCI ~ Foraging Region + DOY + Year
# Creates a table and summarizes the statistics
# Returns the statistics
# Similarly this is hard coded, but could be easily replaced with your covariates that fall in the same places 
extract_etasq_func <- function(list){
  etasq_table <- list[[3]]
  etasq_df <- as.data.frame(etasq_table)
  foraging_region_etasq <- etasq_df$eta.sq[1]
  foraging_region_etasqpart <- etasq_df$eta.sq.part[1]
  DOY_etasq <- etasq_df$eta.sq[2]
  DOY_etasqpart <- etasq_df$eta.sq.part[2]
  Year_etasq <- etasq_df$eta.sq[3]
  Year_etasqpart <- etasq_df$eta.sq.part[3]
  etasq_results_t <- as.data.frame(rbind(foraging_region_etasq, foraging_region_etasqpart, DOY_etasq, DOY_etasqpart, Year_etasq, Year_etasqpart))
  etasq_results <- data.table::transpose(etasq_results_t)
  colnames(etasq_results) <- rownames(etasq_results_t)
  
  return(etasq_results)
  
}

# This function looks at the relationship between DOY and BCI within each region
# First filters by region
# Then creates looks at the linear relationship between BCI and DOY for each of those regions
# Extracts the relevant statistics from each of the linear models for each region
# Pulls all of the relevant statistics into a dataframe 
doy_bci_seasonal <- function(df){
  temp <- df
  
  # Filters by region 
  temp_nyb <- temp %>% filter(Foraging_Region == "New York Bight")
  temp_ice <-  temp %>% filter(Foraging_Region == "Iceland")
  temp_gom <- temp %>% filter(Foraging_Region == "Gulf of Maine")
  temp_gl <-  temp %>% filter(Foraging_Region == "Greenland")
  
  # Calculates BCI ~ DOY for each region and all regions pooled together 
  sum_nyb <- summary(lm(temp_nyb$BCI ~ temp_nyb$DOY))
  sum_ice <- summary(lm(temp_ice$BCI ~ temp_ice$DOY))
  sum_all <- summary(lm(temp$BCI ~ temp$DOY))
  sum_gom <- summary(lm(temp_gom$BCI ~ temp_gom$DOY))
  sum_gl <- summary(lm(temp_gl$BCI ~ temp_gl$DOY))
  
  LM_BCIDOY <- lm(temp_ice$BCI ~ temp_ice$DOY) 
  
  # Pulls out the relevant statistics from the summary(lm) for each region
  nyb_intercept_estimate <- sum_nyb$coefficients[1,1]
  nyb_DOY_estimate <- sum_nyb$coefficients[2,1]
  nyb_intercept_stderror <- sum_nyb$coefficients[2,1]
  nyb_DOY_stderror <- sum_nyb$coefficients[2,2]
  nyb_intercept_tvalue <- sum_nyb$coefficients[1,3]
  nyb_DOY_tvalue <- sum_nyb$coefficients[2,3]
  nyb_intercept_pvalue <- sum_nyb$coefficients[1,4]
  nyb_DOY_pvalue <- sum_nyb$coefficients[2,4]
  
  all_intercept_estimate <- sum_all$coefficients[1,1]
  all_DOY_estimate <- sum_all$coefficients[2,1]
  all_intercept_stderror <- sum_all$coefficients[2,1]
  all_DOY_stderror <- sum_all$coefficients[2,2]
  all_intercept_tvalue <- sum_all$coefficients[1,3]
  all_DOY_tvalue <- sum_all$coefficients[2,3]
  all_intercept_pvalue <- sum_all$coefficients[1,4]
  all_DOY_pvalue <- sum_all$coefficients[2,4]
  
  ice_intercept_estimate <- sum_ice$coefficients[1,1]
  ice_DOY_estimate <- sum_ice$coefficients[2,1]
  ice_intercept_stderror <- sum_ice$coefficients[2,1]
  ice_DOY_stderror <- sum_ice$coefficients[2,2]
  ice_intercept_tvalue <- sum_ice$coefficients[1,3]
  ice_DOY_tvalue <- sum_ice$coefficients[2,3]
  ice_intercept_pvalue <- sum_ice$coefficients[1,4]
  ice_DOY_pvalue <- sum_ice$coefficients[2,4]
  
  ice_intercept <- as.numeric(LM_BCIDOY$coefficients[1])
  ice_slope <- as.numeric(LM_BCIDOY$coefficients[2])
  ice_rsq <- sum_ice$r.squared
  ice_adj_rsq <- sum_ice$adj.r.squared
  
  gom_intercept_estimate <- sum_gom$coefficients[1,1]
  gom_DOY_estimate <- sum_gom$coefficients[2,1]
  gom_intercept_stderror <- sum_gom$coefficients[2,1]
  gom_DOY_stderror <- sum_gom$coefficients[2,2]
  gom_intercept_tvalue <- sum_gom$coefficients[1,3]
  gom_DOY_tvalue <- sum_gom$coefficients[2,3]
  gom_intercept_pvalue <- sum_gom$coefficients[1,4]
  gom_DOY_pvalue <- sum_gom$coefficients[2,4]
  
  gl_intercept_estimate <- sum_gl$coefficients[1,1]
  gl_DOY_estimate <- sum_gl$coefficients[2,1]
  gl_intercept_stderror <- sum_gl$coefficients[2,1]
  gl_DOY_stderror <- sum_gl$coefficients[2,2]
  gl_intercept_tvalue <- sum_gl$coefficients[1,3]
  gl_DOY_tvalue <- sum_gl$coefficients[2,3]
  gl_intercept_pvalue <- sum_gl$coefficients[1,4]
  gl_DOY_pvalue <- sum_gl$coefficients[2,4]
  
  # Binds all of the extracted values into a dataframe
  seasonal_results_t <- as.data.frame(rbind(nyb_intercept_estimate, nyb_DOY_estimate, nyb_intercept_stderror, nyb_DOY_stderror,
                                            nyb_intercept_tvalue, nyb_DOY_tvalue, nyb_intercept_pvalue, nyb_DOY_pvalue,
                                            all_intercept_estimate, all_DOY_estimate, all_intercept_stderror, all_DOY_stderror,
                                            all_intercept_tvalue, all_DOY_tvalue, all_intercept_pvalue, all_DOY_pvalue,
                                            gl_intercept_estimate, gl_DOY_estimate, gl_intercept_stderror, gl_DOY_stderror,
                                            gl_intercept_tvalue, gl_DOY_tvalue, gl_intercept_pvalue, gl_DOY_pvalue,
                                            gom_intercept_estimate, gom_DOY_estimate, gom_intercept_stderror, gom_DOY_stderror,
                                            gom_intercept_tvalue, gom_DOY_tvalue, gom_intercept_pvalue, gom_DOY_pvalue,
                                            ice_intercept_estimate, ice_DOY_estimate, ice_intercept_stderror, ice_DOY_stderror,
                                            ice_intercept_tvalue, ice_DOY_tvalue, ice_intercept_pvalue, ice_DOY_pvalue,
                                            ice_intercept, ice_slope, ice_rsq, ice_adj_rsq))
  seasonal_results <- data.table::transpose(seasonal_results_t)
  colnames(seasonal_results) <- rownames(seasonal_results_t)
  
  return(seasonal_results)
  
  
}

# This function looks at the relationship between BCI and TL
# First calculates the linear model
# Then pulls out all of the relevant statistics from the linear model
# Pulls all of those statistics into a dataframe 
BCI_TL <- function(df){
  
  temp <- df
  sum_all <- summary(lm(temp$BCI ~ temp$Total_Length))
  
  all_intercept_estimate <- sum_all$coefficients[1,1]
  all_TL_estimate <- sum_all$coefficients[2,1]
  all_intercept_stderror <- sum_all$coefficients[2,1]
  all_TL_stderror <- sum_all$coefficients[2,2]
  all_intercept_tvalue <- sum_all$coefficients[1,3]
  all_TL_tvalue <- sum_all$coefficients[2,3]
  all_intercept_pvalue <- sum_all$coefficients[1,4]
  all_TL_pvalue <- sum_all$coefficients[2,4]
  
  results_t <- as.data.frame(rbind(all_intercept_estimate, all_TL_estimate, all_intercept_stderror, all_TL_stderror,
                                   all_intercept_tvalue, all_TL_tvalue, all_intercept_pvalue, all_TL_pvalue))
  results <- data.table::transpose(results_t)
  colnames(results) <- rownames(results_t)
  
  return(results)
  
}

# Run Analysis On Each Bootstrapped Population ------------------------------------------------------------

# So now that the functions have been set up, we run them on all 1000 populations using lapply

# The distinct part of this section of the code is that after we analyze 1000 populations
# What do we do with the statistics? 
# Well we summarize those relevant statistics by pulling out the min, max, mean, confidence interval, and standard deviation
# From here it's up to you to decide which of those min/max/mean/CI/SD statistics are meaningful
# Aka is the standard deviation of all 1000 p-values important? Absolutely not. A standard deviation of a p-value isn't meaningful at all
# But a minimum p-value of 0.10 means that all 1000 p-values are > 0.05, and THAT is meaningful
# Similarly, a maximum p-value of 0.04 means that all 1000 p-values are < 0.05, which is also meaningful
# For the sake of handling so much data at once from every population, and for so many different statistics
# EVERY statistic has min/max/mean/CI/SD calculated in the code but from there you have to interpret what is meaningful
# to the context of 1000 populations and what is not
# For any statistics that are in the middle (let's say the minimum p-value is 0.04, and the max is 0.60, 
# it's important to dig deeper into the direct function output and see how many rows/populations fall above
# or below that 0.05 threshold. If 100 rows/populations fall below the threshold, that's 10% of your data)

# First ensure columns are in numeric and then calculate BCI 
pop_boot_all <- lapply(pop_boot_list, convert_to_numeric, columns = c("Total_Length", "Body_Volume", "DOY"))
pop_boot_BCI <- lapply(pop_boot_all, BCI_func)

# Calculate age class t-test and bind result rows 
pop_boot_age_ttest <- lapply(pop_boot_BCI, ageclass_ttest_func)
pop_boot_age_ttest_results <- bind_rows(pop_boot_age_ttest)  

# Run age model ancova, extract relevant statistics, bind result rows
pop_boot_ancova_agemodel <- lapply(pop_boot_BCI, ancova_agemodel_func) 
pop_boot_ancova_agemodel_results <- lapply(pop_boot_ancova_agemodel, extract_ancova_agemodel_func) 
pop_boot_ancova_agemodel_results <- bind_rows(pop_boot_ancova_agemodel_results)

# Run ANCOVA (BCI ~ Foraging Region + Year + DOY), extract ANCOVA stats, bind rows
# Extract the Tukey stats and bind the rows
# Extract the EtaSq stats and bind the rows 
pop_boot_ancova <- lapply(pop_boot_BCI, ancova_func)
pop_boot_ancova_results <- lapply(pop_boot_ancova, extract_ancova_func)
pop_boot_ancova_results <- bind_rows(pop_boot_ancova_results)
pop_boot_tukey_results <- lapply(pop_boot_ancova, extract_tukey_func)
pop_boot_tukey_results <- bind_rows(pop_boot_tukey_results)
pop_boot_etasq_results <- lapply(pop_boot_ancova, extract_etasq_func)
pop_boot_etasq_results <- bind_rows(pop_boot_etasq_results)

# Here comes the summarizing statistics for each of the 1000 analyses
# Once again, I'm calculating the mean, CI, SD, Min and Max for all of these stats but that does not
# mean that they are meaningful summaries of the statistics for each row!
# I calculate all of them so that it's easier to bulk-calculate and then
# pick out the column and row that makes sense to interpret than code for each and every specific combo

# this code works by calculating the summary statistics for the 1000 rows for each col, which in itself
# represents a statistic or model output that was extracted and bound in a big dataframe "_results"
pop_boot_ttest_agemodel_stats <- list()
for(i in 1:(ncol(pop_boot_age_ttest_results))){
  Mean <- mean(pop_boot_age_ttest_results[,i], na.rm = TRUE)
  Quant97 <- quantile(pop_boot_age_ttest_results[,i], 0.975, na.rm = TRUE)
  Quant2 <-quantile(pop_boot_age_ttest_results[,i], 0.025, na.rm = TRUE)
  SD <- sd(pop_boot_age_ttest_results[,i], na.rm = TRUE)
  Min <- min(pop_boot_age_ttest_results[,i], na.rm = TRUE)
  Max <- max(pop_boot_age_ttest_results[,i], na.rm = TRUE)
  stats <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats) <- colnames(pop_boot_age_ttest_results[i])
  pop_boot_ttest_agemodel_stats[[i]] <- stats
}  

# Bind the summary stats for the t test age models 
ttest_agemodel_stats_df <- do.call(cbind, pop_boot_ttest_agemodel_stats)

# Calculates the summary stats for the Ancova Age Model results
pop_boot_ancova_agemodel_stats <- list()
for(i in 1:(ncol(pop_boot_ancova_agemodel_results))){
  Mean <- mean(pop_boot_ancova_agemodel_results[,i], na.rm = TRUE)
  Quant97 <- quantile(pop_boot_ancova_agemodel_results[,i], 0.975, na.rm = TRUE)
  Quant2 <-quantile(pop_boot_ancova_agemodel_results[,i], 0.025, na.rm = TRUE)
  SD <- sd(pop_boot_ancova_agemodel_results[,i], na.rm = TRUE)
  Min <- min(pop_boot_ancova_agemodel_results[,i], na.rm = TRUE)
  Max <- max(pop_boot_ancova_agemodel_results[,i], na.rm = TRUE)
  stats <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats) <- colnames(pop_boot_ancova_agemodel_results[i])
  pop_boot_ancova_agemodel_stats[[i]] <- stats
}

# Binds the summary stats for the Ancova Age Model results 
ancova_agemodel_stats_df <- do.call(cbind, pop_boot_ancova_agemodel_stats)

# Calculates the summary stats for the ANCOVA results
pop_boot_ancova_stats <- list()
for(i in 1:(ncol(pop_boot_ancova_results))){
  Mean <- mean(pop_boot_ancova_results[,i])
  Quant97 <- quantile(pop_boot_ancova_results[,i], 0.975)
  Quant2 <-quantile(pop_boot_ancova_results[,i], 0.025)
  SD <- sd(pop_boot_ancova_results[,i])
  Min <- min(pop_boot_ancova_results[,i], na.rm = TRUE)
  Max <- max(pop_boot_ancova_results[,i], na.rm = TRUE)
  stats <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats) <- colnames(pop_boot_ancova_results[i])
  pop_boot_ancova_stats[[i]] <- stats
}

# Binds the summary stats for the ANCOVA results 
ancova_stats_df <- do.call(cbind, pop_boot_ancova_stats)

# Calculates the summary stats for the Tukey results
pop_boot_tukey_stats <- list()
for(i in 1:(ncol(pop_boot_tukey_results))){
  Mean <- mean(pop_boot_tukey_results[,i])
  Quant97 <- quantile(pop_boot_tukey_results[,i], 0.975)
  Quant2 <-quantile(pop_boot_tukey_results[,i], 0.025)
  SD <- sd(pop_boot_tukey_results[,i])
  Min <- min(pop_boot_tukey_results[,i], na.rm = TRUE)
  Max <- max(pop_boot_tukey_results[,i], na.rm = TRUE)
  stats <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats) <- colnames(pop_boot_tukey_results[i])
  pop_boot_tukey_stats[[i]] <- stats
}

# Binds the summary stats for the Tukey results
tukey_stats_df <- do.call(cbind, pop_boot_tukey_stats)

# Calculates the summary stats for the Eta Square Results
pop_boot_etasq_stats <- list()
for(i in 1:(ncol(pop_boot_etasq_results))){
  Mean <- mean(pop_boot_etasq_results[,i])
  Quant97 <- quantile(pop_boot_etasq_results[,i], 0.975)
  Quant2 <-quantile(pop_boot_etasq_results[,i], 0.025)
  SD <- sd(pop_boot_etasq_results[,i])
  Min <- min(pop_boot_etasq_results[,i], na.rm = TRUE)
  Max <- max(pop_boot_etasq_results[,i], na.rm = TRUE)
  stats <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats) <- colnames(pop_boot_etasq_results[i])
  pop_boot_etasq_stats[[i]] <- stats
}

# Binds the summary stats for the Eta Square stats 
etasq_stats_df <- do.call(cbind, pop_boot_etasq_stats)


# Calculating the relationship between BCI ~ DOY for each region and all together
pop_boot_seasonalBCI_list <- lapply(pop_boot_BCI, doy_bci_seasonal)
pop_boot_seasonalBCI <- bind_rows(pop_boot_seasonalBCI_list)

# Summarizing the results of the linear models beetween BCI and DOY for each population 
pop_boot_seasonal_BCIDOY <- list()
for(i in 1:(ncol(pop_boot_seasonalBCI))){
  Mean <- mean(pop_boot_seasonalBCI[,i])
  Quant97 <- quantile(pop_boot_seasonalBCI[,i], 0.975)
  Quant2 <-quantile(pop_boot_seasonalBCI[,i], 0.025)
  SD <- sd(pop_boot_seasonalBCI[,i])
  Min <- min(pop_boot_seasonalBCI[,i])
  Max <- max(pop_boot_seasonalBCI[,i])
  stats_ASO <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats_ASO) <- colnames(pop_boot_seasonalBCI[i])
  pop_boot_seasonal_BCIDOY[[i]] <- stats_ASO
}

# Binding the summary stats for the linear models between DOY ~ BCI for each foraging region and all together
seasonalBCIDOY_results <- do.call(cbind, pop_boot_seasonal_BCIDOY)

# Calculating the linear relationship between BCI ~ TL for all populations
BCI_tl <- lapply(pop_boot_BCI, BCI_TL)
pop_boot_BCI_TL <- bind_rows(BCI_tl)

# Calculating the summary stats for BCI ~ TL linear relationship for all populations  
pop_boot_BCITL <- list()
for(i in 1:(ncol(pop_boot_BCI_TL))){
  Mean <- mean(pop_boot_BCI_TL[,i])
  Quant97 <- quantile(pop_boot_BCI_TL[,i], 0.975)
  Quant2 <-quantile(pop_boot_BCI_TL[,i], 0.025)
  SD <- sd(pop_boot_BCI_TL[,i])
  Min <- min(pop_boot_BCI_TL[,i])
  Max <- max(pop_boot_BCI_TL[,i])
  stats_ASO <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats_ASO) <- colnames(pop_boot_BCI_TL[i])
  pop_boot_BCITL[[i]] <- stats_ASO
}

# Binding the summary stats for BCI ~ TL linear relationship for all populations 
BCI_TL_results <- do.call(cbind, pop_boot_BCITL)


# Run Analysis On Each Bootstrapped Population AUGUST SEPT ONLY  ------------------------------------------------------------

# This section is to run analysis on each bootstrapped population but only for data that was collected
# for August and September
# This is to look at the same analysis when only looking at specific months in which data was collected
# Within the Gulf of Maine, Iceland, and New York Bight

# Because of this, it's essentially a new analysis! and that means new functions to use lapply with
# For notation, anything that ends with "ASO" denotes that those functions are for the August September Only analyses
# The functions are not functionally (haha) in most cases, but instead the adjustments are usually just
# adjusting from 4 to 3 regions (no August or Sept data from Greenland)
# which means that the alphabetical order is also adjusted 

# Filters to data only collected between August and September 
augsept_only <- function(df){
  temp <- df
  temp <- temp %>% filter(DOY > 213) %>% filter(DOY < 273) # august 1 to september 30
  return(temp)
}

# Extracts the Tukey output from only the August and September data
# Binds them into a dataframe 
extract_tukey_ASO_func <- function(list){
  tukey_table <- list[[2]]
  iceland_gom_estimate <- as.numeric(tukey_table$test$coefficients[1])
  iceland_gom_stderror <- as.numeric(tukey_table$test$sigma[1])
  iceland_gom_tvalue <- as.numeric(tukey_table$test$tstat[1])
  iceland_gom_pvalue <- as.numeric(tukey_table$test$pvalues[1])
  nyb_gom_estimate <- as.numeric(tukey_table$test$coefficients[2])
  nyb_gom_stderror <- as.numeric(tukey_table$test$sigma[2])
  nyb_gom_tvalue <- as.numeric(tukey_table$test$tstat[2])
  nyb_gom_pvalue <- as.numeric(tukey_table$test$pvalues[2])
  nyb_iceland_estimate <- as.numeric(tukey_table$test$coefficients[3])
  nyb_iceland_stderror <- as.numeric(tukey_table$test$sigma[3])
  nyb_iceland_tvalue <- as.numeric(tukey_table$test$tstat[3])
  nyb_iceland_pvalue <- as.numeric(tukey_table$test$pvalues[3])
  tukey_results_t <- as.data.frame(rbind(iceland_gom_estimate,iceland_gom_stderror, iceland_gom_tvalue,iceland_gom_pvalue,
                                         nyb_gom_estimate,nyb_gom_stderror, nyb_gom_tvalue,nyb_gom_pvalue,
                                         nyb_iceland_estimate,nyb_iceland_stderror,nyb_iceland_tvalue,nyb_iceland_pvalue))
  tukey_results <- data.table::transpose(tukey_results_t)
  colnames(tukey_results) <- rownames(tukey_results_t)
  
  return(tukey_results)
  
}

# Ensure that TL, BV and DOY are numeric 
pop_boot_all <- lapply(pop_boot_list, convert_to_numeric, columns = c("Total_Length", "Body_Volume", "DOY"))

# Filter all data to only include August and September 
pop_boot_all_ASO <- lapply(pop_boot_all, augsept_only)

# Calculate BCI for each individual in each population that has been filtered to only August and Sept data
pop_boot_BCI_ASO <- lapply(pop_boot_all_ASO, BCI_func)

# Run the ancova, extract results, and bind
# Run the Tukey, extract model outputs, and bind
# Run the EtaSquare, extract outputs, and bind 
pop_boot_ancova_ASO <- lapply(pop_boot_BCI_ASO, ancova_func)
pop_boot_ancova_results_ASO <- lapply(pop_boot_ancova_ASO, extract_ancova_func)
pop_boot_ancova_results_ASO <- bind_rows(pop_boot_ancova_results_ASO)
pop_boot_tukey_results_ASO <- lapply(pop_boot_ancova_ASO, extract_tukey_ASO_func)
pop_boot_tukey_results_ASO <- bind_rows(pop_boot_tukey_results_ASO)
pop_boot_etasq_results_ASO <- lapply(pop_boot_ancova_ASO, extract_etasq_func)
pop_boot_etasq_results_ASO <- bind_rows(pop_boot_etasq_results_ASO)

# Calculate the summary statistics for the Ancova stat outputs for August/September data only 
pop_boot_ancova_stats_ASO <- list()
for(i in 1:(ncol(pop_boot_ancova_results_ASO))){
  Mean <- mean(pop_boot_ancova_results_ASO[,i])
  Quant97 <- quantile(pop_boot_ancova_results_ASO[,i], 0.975)
  Quant2 <-quantile(pop_boot_ancova_results_ASO[,i], 0.025)
  SD <- sd(pop_boot_ancova_results_ASO[,i])
  Min <- min(pop_boot_ancova_results_ASO[,i])
  Max <- max(pop_boot_ancova_results_ASO[,i])
  stats_ASO <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats_ASO) <- colnames(pop_boot_ancova_results_ASO[i])
  pop_boot_ancova_stats_ASO[[i]] <- stats_ASO
}

# Bind the summary stats for the Ancova stats for August/September data only 
ancova_stats_df_ASO <- do.call(cbind, pop_boot_ancova_stats_ASO)

# Calculate the summary statistics for the Tukey stat outputs for August/September data only 
pop_boot_tukey_stats_ASO <- list()
for(i in 1:(ncol(pop_boot_tukey_results_ASO))){
  Mean <- mean(pop_boot_tukey_results_ASO[,i])
  Quant97 <- quantile(pop_boot_tukey_results_ASO[,i], 0.975)
  Quant2 <-quantile(pop_boot_tukey_results_ASO[,i], 0.025)
  SD <- sd(pop_boot_tukey_results_ASO[,i])
  Min <- min(pop_boot_tukey_results_ASO[,i])
  Max <- max(pop_boot_tukey_results_ASO[,i])
  stats_ASO <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats_ASO) <- colnames(pop_boot_tukey_results_ASO[i])
  pop_boot_tukey_stats_ASO[[i]] <- stats_ASO
}

# Bind the summary stats for the Tukey stats for the August / September data only 
tukey_stats_df_ASO <- do.call(cbind, pop_boot_tukey_stats_ASO)

# Calculate the summary statistics for the Eta Square stat outputs for August/September data only 
pop_boot_etasq_stats_ASO <- list()
for(i in 1:(ncol(pop_boot_etasq_results_ASO))){
  Mean <- mean(pop_boot_etasq_results_ASO[,i])
  Quant97 <- quantile(pop_boot_etasq_results_ASO[,i], 0.975)
  Quant2 <-quantile(pop_boot_etasq_results_ASO[,i], 0.025)
  SD <- sd(pop_boot_etasq_results_ASO[,i])
  Min <- min(pop_boot_etasq_results_ASO[,i])
  Max <- max(pop_boot_etasq_results_ASO[,i])
  stats_ASO <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats_ASO) <- colnames(pop_boot_etasq_results_ASO[i])
  pop_boot_etasq_stats_ASO[[i]] <- stats_ASO
}

# Bind the summary stats for the EtaSq for the August / September data only 
etasq_stats_df_ASO <- do.call(cbind, pop_boot_etasq_stats_ASO)





# SMA  --------------------------------------------------------------------

# In this section I calculate an SMA for each of the populations
# The backbone of this section is the smatr package

# Generate the SMA between Body Volume ~ Total Length and Foraging Region
generate_sma_func <- function(df){
  df$Foraging_Region <- gsub("New York Bight", "New_York_Bight", df$Foraging_Region) # important because the smatr package does not do spaces or ` ` for columns 
  df$Foraging_Region <- gsub("Gulf of Maine", "Gulf_of_Maine", df$Foraging_Region) # important because the smatr package does not do spaces or ` ` for columns 
  df$Foraging_Region <- as.factor(df$Foraging_Region) # also important to make it a factor 
  df <- as.data.frame(df)
  sma_result <- smatr::sma(Body_Volume ~ Total_Length + Foraging_Region, log = "xy", data = df, multcomp = TRUE, multcompmethod = "adjusted") 
  return(sma_result)
}

# Generate SMA 
sma_list_no_filter <- lapply(pop_boot_all, generate_sma_func)

# Filter By SMA that have comparable slopes 
sma_list <- c()  # Initialize an empty list
for (i in 1:length(sma_list_no_filter)) {
  slope <- sma_list_no_filter[[i]]$commoncoef$p
  if (slope > 0.05) {
    sma_list[[i]] <- sma_list_no_filter[[i]]  # Add value to the list
  }
}

sma_list <- sma_list[!sapply(sma_list, is.null)]

# The SMA Results are summarized sligthly different here, as a loop for each column
# They are summarized in this huge data table here

sma_RS_results_list <- list()
for(i in 1:length(sma_list)){
  SMA_Index <- i
  H0_Equal_Slope_P <- sma_list[[i]]$commoncoef$p
  H0_Equal_Elev_P <- sma_list[[i]]$gtr$p[1]
  Wald_Statistic_Elev <- sma_list[[i]]$gtr$stat[1]
  Greenland_R2 <- sma_list[[i]]$groupsummary$r2[1]
  Greenland_R2_P <- sma_list[[i]]$groupsummary$pval[1]
  Greenland_Slope <- sma_list[[i]]$groupsummary$Slope[1]
  Greenland_Elevation <- sma_list[[i]]$gtr$as[1]
  Greenland_Plot_Elevation <- sma_list[[i]]$coef$Greenland$`coef(SMA)`[1]
  Greenland_Plot_Slope <- sma_list[[i]]$coef$Greenland$`coef(SMA)`[2]
  Greenland_Plot_Int <- sma_list[[i]]$groupsummary$Int[1]
  GOM_R2 <- sma_list[[i]]$groupsummary$r2[2]
  GOM_R2_P <- sma_list[[i]]$groupsummary$pval[2]
  GOM_Slope <- sma_list[[i]]$groupsummary$Slope[2]
  GOM_Elevation <- sma_list[[i]]$gtr$as[2]
  GOM_Plot_Elevation <- sma_list[[i]]$coef$Gulf_of_Maine$`coef(SMA)`[1]
  GOM_Plot_Slope <- sma_list[[i]]$coef$Gulf_of_Maine$`coef(SMA)`[2]
  GOM_Plot_Int <- sma_list[[i]]$groupsummary$Int[2]
  Iceland_R2 <- sma_list[[i]]$groupsummary$r2[3]
  Iceland_R2_P <- sma_list[[i]]$groupsummary$pval[3]
  Iceland_Slope <- sma_list[[i]]$groupsummary$Slope[3]
  Iceland_Elevation <- sma_list[[i]]$gtr$as[3]
  Iceland_Plot_Elevation <- sma_list[[i]]$coef$Iceland$`coef(SMA)`[1]
  Iceland_Plot_Slope <- sma_list[[i]]$coef$Iceland$`coef(SMA)`[2]
  Iceland_Plot_Int <- sma_list[[i]]$groupsummary$Int[3]
  NYB_R2 <- sma_list[[i]]$groupsummary$r2[4]
  NYB_R2_P <- sma_list[[i]]$groupsummary$pval[4]
  NYB_Slope <- sma_list[[i]]$groupsummary$Slope[4]
  NYB_Elevation <- sma_list[[i]]$gtr$as[4]
  NYB_Plot_Elevation <- sma_list[[i]]$coef$New_York_Bight$`coef(SMA)`[1]
  NYB_Plot_Slope <- sma_list[[i]]$coef$New_York_Bight$`coef(SMA)`[2]
  NYB_Plot_Int <- sma_list[[i]]$groupsummary$Int[4]
  
  GL_GOM_Pval <- sma_list[[i]]$multcompresult$Pval[1]
  GL_GOM_TestStat <- sma_list[[i]]$multcompresult$TestStat[1]
  GL_GOM_Elev1 <- sma_list[[i]]$multcompresult$Elev1[1]
  GL_GOM_Elev2 <- sma_list[[i]]$multcompresult$Elev2[1]
  
  GL_ICE_Pval <- sma_list[[i]]$multcompresult$Pval[2]
  GL_ICE_TestStat <- sma_list[[i]]$multcompresult$TestStat[2]
  GL_ICE_Elev1 <- sma_list[[i]]$multcompresult$Elev1[2]
  GL_ICE_Elev2 <- sma_list[[i]]$multcompresult$Elev2[2]
  
  GL_NYB_Pval <- sma_list[[i]]$multcompresult$Pval[3]
  GL_NYB_TestStat <- sma_list[[i]]$multcompresult$TestStat[3]
  GL_NYB_Elev1 <- sma_list[[i]]$multcompresult$Elev1[3]
  GL_NYB_Elev2 <- sma_list[[i]]$multcompresult$Elev2[3]
  
  ICE_GOM_Pval <- sma_list[[i]]$multcompresult$Pval[4]
  ICE_GOM_TestStat <- sma_list[[i]]$multcompresult$TestStat[4]
  ICE_GOM_Elev1 <- sma_list[[i]]$multcompresult$Elev1[4]
  ICE_GOM_Elev2 <- sma_list[[i]]$multcompresult$Elev2[4]
  
  NYB_GOM_Pval <- sma_list[[i]]$multcompresult$Pval[5]
  NYB_GOM_TestStat <- sma_list[[i]]$multcompresult$TestStat[5]
  NYB_GOM_Elev1 <- sma_list[[i]]$multcompresult$Elev1[5]
  NYB_GOM_Elev2 <- sma_list[[i]]$multcompresult$Elev2[5]
  
  NYB_ICE_Pval <- sma_list[[i]]$multcompresult$Pval[6]
  NYB_ICE_TestStat <- sma_list[[i]]$multcompresult$TestStat[6]
  NYB_ICE_Elev1 <- sma_list[[i]]$multcompresult$Elev1[6]
  NYB_ICE_Elev2 <- sma_list[[i]]$multcompresult$Elev2[6]
  
  
  temp_sma_results_df <- as.data.frame(cbind(SMA_Index, H0_Equal_Slope_P, H0_Equal_Elev_P, Wald_Statistic_Elev, 
                                             Greenland_R2, Greenland_R2_P, Greenland_Slope, Greenland_Elevation, Greenland_Plot_Elevation, Greenland_Plot_Slope, Greenland_Plot_Int,
                                             GOM_R2, GOM_R2_P, GOM_Slope, GOM_Elevation, GOM_Plot_Elevation, GOM_Plot_Slope, GOM_Plot_Int, 
                                             Iceland_R2, Iceland_R2_P, Iceland_Slope, Iceland_Elevation, Iceland_Plot_Elevation, Iceland_Plot_Slope, Iceland_Plot_Int,
                                             NYB_R2, NYB_R2_P, NYB_Slope, NYB_Elevation, NYB_Plot_Elevation, NYB_Plot_Slope, NYB_Plot_Int,
                                             GL_GOM_Elev1, GL_GOM_Elev2, GL_GOM_Pval, GL_GOM_TestStat, GL_ICE_Elev1, GL_ICE_Elev2, GL_ICE_Pval,     
                                             GL_ICE_TestStat, GL_NYB_Elev1, GL_NYB_Elev2, GL_NYB_Pval, GL_NYB_TestStat, ICE_GOM_Elev1, ICE_GOM_Elev2,   
                                             ICE_GOM_Pval, ICE_GOM_TestStat, NYB_GOM_Elev1, NYB_GOM_Elev2, NYB_GOM_Pval, NYB_GOM_TestStat, NYB_ICE_Elev1,   
                                             NYB_ICE_Elev2, NYB_ICE_Pval, NYB_ICE_TestStat
  ))
  sma_RS_results_list[[i]] <- temp_sma_results_df
}


# These are all of the results of all 1000 of the SMAs
# They are bound here 
sma_RS_results_DF <- do.call(rbind, c(sma_RS_results_list, make.row.names = FALSE))

# Collate all of the raw data points used in the SMA so that you can plot it later
# This is important because the SMA calculates the log10 of all of the points, otherwise you'd have to do it
# in another step
# Also the smatr package doesn't have compatability with ggplot, so this is important to do it manually
sma_data_list <- list()
for(i in 1:length(sma_list)){
  sma_data_frame <- sma_list[[i]]$data
  sma_data_list[[i]] <- sma_data_frame
}

# Collated SMA$Data from all 1000 individuals 
sma_collated_data <- do.call(rbind, c(sma_data_list, make.row.names = FALSE))

# A Loop to Find All of the Summary Stats Related to the Bootstrapped SMAs
# Same thing as above - all of the mean/CI/SD/min/max stats are calculated for all of the SMA stats
# But that doesn't mean that they are meaningful! Interpret carefully which summary stats are important for each SMA stat 

stats_list <- list()

for(i in 1:(ncol(sma_RS_results_DF))){
  Mean <- mean(sma_RS_results_DF[,i])
  Quant97 <- quantile(sma_RS_results_DF[,i], 0.975)
  Quant2 <-quantile(sma_RS_results_DF[,i], 0.025)
  SD <- sd(sma_RS_results_DF[,i])
  Min <- min(sma_RS_results_DF[,i])
  Max <- max(sma_RS_results_DF[,i])
  stats <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats) <- colnames(sma_RS_results_DF[i])
  stats_list[[i]] <- stats
}

# Bind all of the SMA summary stats into a single data frame 
sma_stats_df <- do.call(cbind, stats_list) %>% dplyr::select(-SMA_Index)

# SMA Mass ----------------------------------------------------------------

# Calculating the SMA Mass is nearly identical to the SMA section, only adding in Mass as a multiplier

# The important first step is to generate a density estimate to multiply volume by
# The density estimate here incorporates the confidence interval around the mean density estimate

generate_sma_mass_func <- function(df){
  
  df$r_density <- rnorm(1, mean = 1037.2, sd = 11.75) # Randomly select a density estimate that is centered around the mean and standard deviation from Aoki et al
  
  df$mass <- df$Body_Volume * df$r_density # multiply density and volume to get mass 
  
  df$Foraging_Region <- gsub("New York Bight", "New_York_Bight", df$Foraging_Region) # important because no spaces allowed in SMA columns
  df$Foraging_Region <- gsub("Gulf of Maine", "Gulf_of_Maine", df$Foraging_Region) # important because no spaces allowed in SMA columns 
  df$Foraging_Region <- as.factor(df$Foraging_Region)
  df <- as.data.frame(df)
  sma_result <- smatr::sma(mass ~ Total_Length + Foraging_Region, log = "xy", data = df, multcomp = TRUE, multcompmethod = "adjusted") 
  return(sma_result)
}

# Generate SMA 
sma_list_mass_no_filter <- lapply(pop_boot_all, generate_sma_mass_func)

# Filter any SMA that do not have comparable slopes 
sma_list_mass <- c()  # Initialize an empty list
for (i in 1:length(sma_list_mass_no_filter)) {
  slope <- sma_list_mass_no_filter[[i]]$commoncoef$p
  if (slope > 0.05) {
    sma_list_mass[[i]] <- sma_list_mass_no_filter[[i]]  # Add value to the list
  }
}

sma_list_mass <- sma_list_mass[!sapply(sma_list_mass, is.null)]

# Collate SMA Results 

sma_RS_results_list_mass <- list()
for(i in 1:length(sma_list_mass)){
  SMA_Index <- i
  H0_Equal_Slope_P <- sma_list_mass[[i]]$commoncoef$p
  H0_Equal_Elev_P <- sma_list_mass[[i]]$gtr$p[1]
  Wald_Statistic_Elev <- sma_list_mass[[i]]$gtr$stat[1]
  Greenland_R2 <- sma_list_mass[[i]]$groupsummary$r2[1]
  Greenland_R2_P <- sma_list_mass[[i]]$groupsummary$pval[1]
  Greenland_Slope <- sma_list_mass[[i]]$groupsummary$Slope[1]
  Greenland_Elevation <- sma_list_mass[[i]]$gtr$as[1]
  Greenland_Plot_Elevation <- sma_list_mass[[i]]$coef$Greenland$`coef(SMA)`[1]
  Greenland_Plot_Slope <- sma_list_mass[[i]]$coef$Greenland$`coef(SMA)`[2]
  Greenland_Plot_Int <- sma_list_mass[[i]]$groupsummary$Int[1]
  GOM_R2 <- sma_list_mass[[i]]$groupsummary$r2[2]
  GOM_R2_P <- sma_list_mass[[i]]$groupsummary$pval[2]
  GOM_Slope <- sma_list_mass[[i]]$groupsummary$Slope[2]
  GOM_Elevation <- sma_list_mass[[i]]$gtr$as[2]
  GOM_Plot_Elevation <- sma_list_mass[[i]]$coef$Gulf_of_Maine$`coef(SMA)`[1]
  GOM_Plot_Slope <- sma_list_mass[[i]]$coef$Gulf_of_Maine$`coef(SMA)`[2]
  GOM_Plot_Int <- sma_list_mass[[i]]$groupsummary$Int[2]
  Iceland_R2 <- sma_list_mass[[i]]$groupsummary$r2[3]
  Iceland_R2_P <- sma_list_mass[[i]]$groupsummary$pval[3]
  Iceland_Slope <- sma_list_mass[[i]]$groupsummary$Slope[3]
  Iceland_Elevation <- sma_list_mass[[i]]$gtr$as[3]
  Iceland_Plot_Elevation <- sma_list_mass[[i]]$coef$Iceland$`coef(SMA)`[1]
  Iceland_Plot_Slope <- sma_list_mass[[i]]$coef$Iceland$`coef(SMA)`[2]
  Iceland_Plot_Int <- sma_list_mass[[i]]$groupsummary$Int[3]
  NYB_R2 <- sma_list_mass[[i]]$groupsummary$r2[4]
  NYB_R2_P <- sma_list_mass[[i]]$groupsummary$pval[4]
  NYB_Slope <- sma_list_mass[[i]]$groupsummary$Slope[4]
  NYB_Elevation <- sma_list_mass[[i]]$gtr$as[4]
  NYB_Plot_Elevation <- sma_list_mass[[i]]$coef$New_York_Bight$`coef(SMA)`[1]
  NYB_Plot_Slope <- sma_list_mass[[i]]$coef$New_York_Bight$`coef(SMA)`[2]
  NYB_Plot_Int <- sma_list_mass[[i]]$groupsummary$Int[4]
  
  GL_GOM_Pval <- sma_list_mass[[i]]$multcompresult$Pval[1]
  GL_GOM_TestStat <- sma_list_mass[[i]]$multcompresult$TestStat[1]
  GL_GOM_Elev1 <- sma_list_mass[[i]]$multcompresult$Elev1[1]
  GL_GOM_Elev2 <- sma_list_mass[[i]]$multcompresult$Elev2[1]
  
  GL_ICE_Pval <- sma_list_mass[[i]]$multcompresult$Pval[2]
  GL_ICE_TestStat <- sma_list_mass[[i]]$multcompresult$TestStat[2]
  GL_ICE_Elev1 <- sma_list_mass[[i]]$multcompresult$Elev1[2]
  GL_ICE_Elev2 <- sma_list_mass[[i]]$multcompresult$Elev2[2]
  
  GL_NYB_Pval <- sma_list_mass[[i]]$multcompresult$Pval[3]
  GL_NYB_TestStat <- sma_list_mass[[i]]$multcompresult$TestStat[3]
  GL_NYB_Elev1 <- sma_list_mass[[i]]$multcompresult$Elev1[3]
  GL_NYB_Elev2 <- sma_list_mass[[i]]$multcompresult$Elev2[3]
  
  ICE_GOM_Pval <- sma_list_mass[[i]]$multcompresult$Pval[4]
  ICE_GOM_TestStat <- sma_list_mass[[i]]$multcompresult$TestStat[4]
  ICE_GOM_Elev1 <- sma_list_mass[[i]]$multcompresult$Elev1[4]
  ICE_GOM_Elev2 <- sma_list_mass[[i]]$multcompresult$Elev2[4]
  
  NYB_GOM_Pval <- sma_list_mass[[i]]$multcompresult$Pval[5]
  NYB_GOM_TestStat <- sma_list_mass[[i]]$multcompresult$TestStat[5]
  NYB_GOM_Elev1 <- sma_list_mass[[i]]$multcompresult$Elev1[5]
  NYB_GOM_Elev2 <- sma_list_mass[[i]]$multcompresult$Elev2[5]
  
  NYB_ICE_Pval <- sma_list_mass[[i]]$multcompresult$Pval[6]
  NYB_ICE_TestStat <- sma_list_mass[[i]]$multcompresult$TestStat[6]
  NYB_ICE_Elev1 <- sma_list_mass[[i]]$multcompresult$Elev1[6]
  NYB_ICE_Elev2 <- sma_list_mass[[i]]$multcompresult$Elev2[6]
  
  
  temp_sma_results_df <- as.data.frame(cbind(SMA_Index, H0_Equal_Slope_P, H0_Equal_Elev_P, Wald_Statistic_Elev, 
                                             Greenland_R2, Greenland_R2_P, Greenland_Slope, Greenland_Elevation, Greenland_Plot_Elevation, Greenland_Plot_Slope, Greenland_Plot_Int,
                                             GOM_R2, GOM_R2_P, GOM_Slope, GOM_Elevation, GOM_Plot_Elevation, GOM_Plot_Slope, GOM_Plot_Int, 
                                             Iceland_R2, Iceland_R2_P, Iceland_Slope, Iceland_Elevation, Iceland_Plot_Elevation, Iceland_Plot_Slope, Iceland_Plot_Int,
                                             NYB_R2, NYB_R2_P, NYB_Slope, NYB_Elevation, NYB_Plot_Elevation, NYB_Plot_Slope, NYB_Plot_Int,
                                             GL_GOM_Elev1, GL_GOM_Elev2, GL_GOM_Pval, GL_GOM_TestStat, GL_ICE_Elev1, GL_ICE_Elev2, GL_ICE_Pval,     
                                             GL_ICE_TestStat, GL_NYB_Elev1, GL_NYB_Elev2, GL_NYB_Pval, GL_NYB_TestStat, ICE_GOM_Elev1, ICE_GOM_Elev2,   
                                             ICE_GOM_Pval, ICE_GOM_TestStat, NYB_GOM_Elev1, NYB_GOM_Elev2, NYB_GOM_Pval, NYB_GOM_TestStat, NYB_ICE_Elev1,   
                                             NYB_ICE_Elev2, NYB_ICE_Pval, NYB_ICE_TestStat
  ))
  sma_RS_results_list_mass[[i]] <- temp_sma_results_df
}


# These are all of the results of all 1000 of the SMAs 
# Bind them 
sma_RS_results_DF_mass <- do.call(rbind, c(sma_RS_results_list_mass, make.row.names = FALSE))

# Collate all of the raw data points used in the SMA so that you can plot it later
# Again, this would need to be done anyway because of the log10 calculation
# and because smatr is not compatible with ggplot and i'm not using Base R plotting (no thank you)
sma_data_list_mass <- list()
for(i in 1:length(sma_list_mass)){
  sma_data_frame <- sma_list_mass[[i]]$data
  sma_data_frame <- sma_data_frame %>% dplyr::select(-Foraging_Region)
  sma_data_frame$ID <- seq.int(nrow(sma_data_frame))
  sma_data_list_mass[[i]] <- sma_data_frame
}

# Collated SMA$Data from all 1000 individuals 
sma_collated_data_mass <- do.call(rbind, c(sma_data_list_mass, make.row.names = FALSE))


# Mean of SMA Mass Data

mean_mass_data <- sma_collated_data_mass %>%
  group_by(ID) %>%
  summarise_all(mean)

mean_mass_data$Foraging_Region <- sma_list_mass[[1]]$data$Foraging_Region


# A Loop to Find All of the Stats Related to the Bootstrapped SMAs
# Summary stats once again 

stats_list <- list()

for(i in 1:(ncol(sma_RS_results_DF_mass))){
  Mean <- mean(sma_RS_results_DF_mass[,i])
  Quant97 <- quantile(sma_RS_results_DF_mass[,i], 0.975)
  Quant2 <-quantile(sma_RS_results_DF_mass[,i], 0.025)
  SD <- sd(sma_RS_results_DF_mass[,i])
  Min <- min(sma_RS_results_DF_mass[,i])
  Max <- max(sma_RS_results_DF_mass[,i])
  stats <- as.data.frame(rbind(Mean, Quant97, Quant2, SD, Min, Max))
  colnames(stats) <- colnames(sma_RS_results_DF_mass[i])
  stats_list[[i]] <- stats
}

sma_stats_mass_df <- do.call(cbind, stats_list) %>% dplyr::select(-SMA_Index)


# Plotting BCI Plots (Figure 1, Figure 3a and b) ------------------------------------------------------------

# Find the Mean BCI for each individual whale

megapop_BCI_df <- bind_rows(pop_boot_BCI) # this just makes it one huge data frame

all_Animal_ID <- unique(pop_boot_BCI[[1]]$Animal_ID) # all the unique Animal_IDs 

# this dataframe takes all of the measurements that incorporate error for each unique Animal_ID
# and finds the mean values of TL, BV, BCI and LogTL, LogBV so that it's easier to plot 
mean_animal_df <- data.frame()
for (i in 1:length(all_Animal_ID)){
  animal_BCI_df <- megapop_BCI_df %>% filter(Animal_ID == all_Animal_ID[i])
  animal_BCI_ref_row <- animal_BCI_df[1,]
  animal_BCI_ref_row$Mean_BCI <-  mean(animal_BCI_df$BCI)
  animal_BCI_ref_row$SD_BCI <- sd(animal_BCI_df$BCI)
  animal_BCI_ref_row$Mean_TL <- mean(animal_BCI_df$Total_Length)
  animal_BCI_ref_row$Mean_BV <- mean(animal_BCI_df$Body_Volume)
  animal_BCI_ref_row$LogMean_TL <- log10(animal_BCI_ref_row$Mean_TL)
  animal_BCI_ref_row$LogMean_BV <- log10(animal_BCI_ref_row$Mean_BV)
  mean_animal_df <- as.data.frame(rbind(animal_BCI_ref_row, mean_animal_df))
}

# Violin Plot Using Mean BCI Values for Each Whale 

mean_animal_df$Foraging_Region <- factor(mean_animal_df$Foraging_Region, levels = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland"))

site_colors <- c("New York Bight" = "firebrick3", "Iceland" = "lightskyblue", "Gulf of Maine" = "goldenrod2", "Historical" = "navy", "Greenland" = "forestgreen")

BCI_violin <- ggplot(mean_animal_df, aes(x = Foraging_Region, y = Mean_BCI, fill = Foraging_Region, alpha = 0.75, color = Foraging_Region)) + # all regional data, with foraging region as x groups and BCI on y axis. Colors determined by foraging group
  geom_violin(fill = "transparent", draw_quantiles = 0.5, color = "black") + # formatting for violin plot
  geom_violin() + # telling ggplot to make a violin plot
  geom_dotplot(binaxis = 'y', stackratio = 1, stackgroups = TRUE, stackdir = 'centerwhole', dotsize = 0.3, alpha = 2, binpositions = "bygroup", stroke = NA, method = "histodot") +
  labs(x = "Foraging Region", y = "Body Condition Index") + # labels
  #ggtitle("Body Condition Index by Region") + # title
  scale_fill_manual(values = site_colors) + # legend 
  scale_color_manual(values = site_colors) +  # legend 
  scale_alpha(guide = "none") +  # legend
  guides(fill="none") +
  guides(color = "none") +
  theme(legend.position = "none") + # title formatting
  geom_signif(comparisons = list(c("Gulf of Maine", "New York Bight")), map_signif_level=TRUE, textsize = 5, fontface = "bold", y_position =  0.62, aes(color = "black"), alpha = 0.8) + # making significance stars between NYB and GOM
  geom_signif(comparisons = list(c("Gulf of Maine", "Iceland")), map_signif_level=TRUE, textsize = 5, fontface = "bold", y_position =  0.56, aes(color = "black"), alpha = 0.8) + # making significance stars between GOM and Iceland
  geom_signif(comparisons = list(c("Gulf of Maine", "Greenland")), map_signif_level=TRUE, textsize = 5, fontface = "bold", y_position =  0.49, aes(color = "black"), alpha = 0.8) + # making significance stars between NYB and Greenland
  theme_bw() + # making theme black and white, removing the default grey background 
  theme(axis.text=element_text(size=12), title=element_text(size=12))
BCI_violin # plot the plot 

# TL v BV Plot

bci_equation <- lm(log(mean_animal_df$Mean_BV) ~ log(mean_animal_df$Mean_TL)) # equation to find expected body volume from Christiansen et al (2020)

log_bci_equation <- lm(log10(mean_animal_df$Mean_BV) ~ log10(mean_animal_df$Mean_TL)) # equation to find expected body volume from Christiansen et al (2020)

vol_by_length <- nls(Mean_BV ~ a*(Mean_TL^b),
                     data = mean_animal_df, 
                     start = list(a = 0.5, b = 0.2))

tl_v_bv_plot <- ggplot(mean_animal_df, aes(x = Mean_TL, y = Mean_BV)) + # all historical data with total length as x and body volume as y 
  geom_point(data = mean_animal_df, aes(x = Mean_TL, y = Mean_BV, fill = Foraging_Region, color = Foraging_Region), size = 2, alpha = 0.75) +
  geom_line(aes(y = exp(fitted(bci_equation)))) + # adding in the expected value line
  labs(x = "Total Length (m)", y = bquote('Body Volume '(m^3)), color = "Legend") + # title and labels 
  scale_fill_manual(values = site_colors, name = "Foraging_Region", breaks = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland")) + # legend 
  scale_color_manual(values = site_colors, name = "Foraging_Region", breaks = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland")) + # legend 
  theme_bw() # + # theme changed to black and white to get rid of the default ggplot 

tl_v_bv_plot + theme(legend.position = "none", 
                     legend.text=element_text(size=19), 
                     legend.title = element_text(size=19),
                     plot.subtitle=element_text(size=22),
                     axis.title = element_text(size=22),
                     axis.text = element_text(size=19),
                     plot.title=element_text(face = "bold", size = 22)) + guides(color = guide_legend(nrow = 3))

log_tl_v_bv_plot <- ggplot(mean_animal_df, aes(x = LogMean_TL, y = LogMean_BV)) + # all historical data with total length as x and body volume as y 
  geom_point(data = mean_animal_df, aes(x = LogMean_TL, y = LogMean_BV, fill = Foraging_Region, color = Foraging_Region), size = 2, alpha = 0.75) +
  geom_line(aes(y = (fitted(log_bci_equation)))) + # adding in the expected value line
  labs(x = "Log Total Length (m)", y = bquote('Log Body Volume '(m^3)), color = "Legend") + # title and labels 
  scale_fill_manual(values = site_colors, name = "Foraging_Region", breaks = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland")) + # legend 
  scale_color_manual(values = site_colors, name = "Foraging_Region", breaks = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland")) + # legend 
  theme_bw() # + # theme changed to black and white to get rid of the default ggplot 

log_tl_v_bv_plot + theme(legend.position = "none", 
                         legend.text=element_text(size=19), 
                         legend.title = element_text(size=19),
                         plot.subtitle=element_text(size=22),
                         axis.title = element_text(size=22),
                         axis.text = element_text(size=19),
                         plot.title=element_text(face = "bold", size = 22)) + guides(color = guide_legend(nrow = 3))  +
  coord_cartesian(xlim = c(0.81, 1.19), ylim = c(0.6, 1.76))


# Plotting The SMA (Figure 3c) !  -------------------------------------------------------------

# Buckle up people - the SMA package, smatr, and ggplot are not compatible
# So what do we do? Use Base R plots? ... absolutely not. 
# Instead I have looked under the hood of how ggsmatr works (a package that attempts to plot smas with ggplot but doesn't work with log log plots)
# And now this convulated plotting code exists
# the same note applies to plotting the SMA for Mass 

# SMA Plot 

colnames_plot <- c("group", "elevation", "slope", "min_x", "max_x", "min_y","max_y", "r2", "pval", "Slope", "Int")

# The following code pulls out the plotting information for each Region mean elevation, and the elevation at 97.5% and 2.5% to form the confidence bands of each

Greenland_Mean <- c("Greenland_Mean", 
                    sma_stats_df$Greenland_Plot_Elevation[1], 
                    sma_stats_df$Greenland_Plot_Slope[1], 
                    min((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                    max((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                    (sma_stats_df$Greenland_Plot_Slope[1]*min((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_df$Greenland_Plot_Elevation[1],
                    (sma_stats_df$Greenland_Plot_Slope[1]*max((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_df$Greenland_Plot_Elevation[1],
                    sma_stats_df$Greenland_R2[1],
                    sma_stats_df$Greenland_R2_P[1],
                    sma_stats_df$Greenland_Slope[1],
                    sma_stats_df$Greenland_Plot_Int[1])

Greenland_97 <- c("Greenland_97", 
                  sma_stats_df$Greenland_Plot_Elevation[2], 
                  sma_stats_df$Greenland_Plot_Slope[2], 
                  min((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                  max((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                  (sma_stats_df$Greenland_Plot_Slope[2]*min((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_df$Greenland_Plot_Elevation[2],
                  (sma_stats_df$Greenland_Plot_Slope[2]*max((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_df$Greenland_Plot_Elevation[2],
                  sma_stats_df$Greenland_R2[2],
                  sma_stats_df$Greenland_R2_P[2],
                  sma_stats_df$Greenland_Slope[2],
                  sma_stats_df$Greenland_Plot_Int[2])

Greenland_2 <- c("Greenland_2", 
                 sma_stats_df$Greenland_Plot_Elevation[3], 
                 sma_stats_df$Greenland_Plot_Slope[3], 
                 min((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                 max((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                 (sma_stats_df$Greenland_Plot_Slope[3]*min((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_df$Greenland_Plot_Elevation[3],
                 (sma_stats_df$Greenland_Plot_Slope[3]*max((sma_collated_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_df$Greenland_Plot_Elevation[3],
                 sma_stats_df$Greenland_R2[3],
                 sma_stats_df$Greenland_R2_P[3],
                 sma_stats_df$Greenland_Slope[3],
                 sma_stats_df$Greenland_Plot_Int[3])

Iceland_Mean <- c("Iceland_Mean", 
                  sma_stats_df$Iceland_Plot_Elevation[1], 
                  sma_stats_df$Iceland_Plot_Slope[1], 
                  min((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                  max((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                  (sma_stats_df$Iceland_Plot_Slope[1]*min((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_df$Iceland_Plot_Elevation[1],
                  (sma_stats_df$Iceland_Plot_Slope[1]*max((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_df$Iceland_Plot_Elevation[1],
                  sma_stats_df$Iceland_R2[1],
                  sma_stats_df$Iceland_R2_P[1],
                  sma_stats_df$Iceland_Slope[1],
                  sma_stats_df$Iceland_Plot_Int[1])

Iceland_97 <- c("Iceland_97", 
                sma_stats_df$Iceland_Plot_Elevation[2], 
                sma_stats_df$Iceland_Plot_Slope[2], 
                min((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                max((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                (sma_stats_df$Iceland_Plot_Slope[2]*min((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_df$Iceland_Plot_Elevation[2],
                (sma_stats_df$Iceland_Plot_Slope[2]*max((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_df$Iceland_Plot_Elevation[2],
                sma_stats_df$Iceland_R2[2],
                sma_stats_df$Iceland_R2_P[2],
                sma_stats_df$Iceland_Slope[2],
                sma_stats_df$Iceland_Plot_Int[2])

Iceland_2 <- c("Iceland_2", 
               sma_stats_df$Iceland_Plot_Elevation[3], 
               sma_stats_df$Iceland_Plot_Slope[3], 
               min((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
               max((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
               (sma_stats_df$Iceland_Plot_Slope[3]*min((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_df$Iceland_Plot_Elevation[3],
               (sma_stats_df$Iceland_Plot_Slope[3]*max((sma_collated_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_df$Iceland_Plot_Elevation[3],
               sma_stats_df$Iceland_R2[3],
               sma_stats_df$Iceland_R2_P[3],
               sma_stats_df$Iceland_Slope[3],
               sma_stats_df$Iceland_Plot_Int[3])

NYB_Mean <- c("NYB_Mean", 
              sma_stats_df$NYB_Plot_Elevation[1], 
              sma_stats_df$NYB_Plot_Slope[1], 
              min((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
              max((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
              (sma_stats_df$NYB_Plot_Slope[1]*min((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_df$NYB_Plot_Elevation[1],
              (sma_stats_df$NYB_Plot_Slope[1]*max((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_df$NYB_Plot_Elevation[1],
              sma_stats_df$NYB_R2[1],
              sma_stats_df$NYB_R2_P[1],
              sma_stats_df$NYB_Slope[1],
              sma_stats_df$NYB_Plot_Int[1])

NYB_97 <- c("NYB_97", 
            sma_stats_df$NYB_Plot_Elevation[2], 
            sma_stats_df$NYB_Plot_Slope[2], 
            min((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
            max((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
            (sma_stats_df$NYB_Plot_Slope[2]*min((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_df$NYB_Plot_Elevation[2],
            (sma_stats_df$NYB_Plot_Slope[2]*max((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_df$NYB_Plot_Elevation[2],
            sma_stats_df$NYB_R2[2],
            sma_stats_df$NYB_R2_P[2],
            sma_stats_df$NYB_Slope[2],
            sma_stats_df$NYB_Plot_Int[2])

NYB_2 <- c("NYB_2", 
           sma_stats_df$NYB_Plot_Elevation[3], 
           sma_stats_df$NYB_Plot_Slope[3], 
           min((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
           max((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
           (sma_stats_df$NYB_Plot_Slope[3]*min((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_df$NYB_Plot_Elevation[3],
           (sma_stats_df$NYB_Plot_Slope[3]*max((sma_collated_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_df$NYB_Plot_Elevation[3],
           sma_stats_df$NYB_R2[3],
           sma_stats_df$NYB_R2_P[3],
           sma_stats_df$NYB_Slope[3],
           sma_stats_df$NYB_Plot_Int[3])

GOM_Mean <- c("GOM_Mean", 
              sma_stats_df$GOM_Plot_Elevation[1], 
              sma_stats_df$GOM_Plot_Slope[1], 
              min((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
              max((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
              (sma_stats_df$GOM_Plot_Slope[1]*min((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_df$GOM_Plot_Elevation[1],
              (sma_stats_df$GOM_Plot_Slope[1]*max((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_df$GOM_Plot_Elevation[1],
              sma_stats_df$GOM_R2[1],
              sma_stats_df$GOM_R2_P[1],
              sma_stats_df$GOM_Slope[1],
              sma_stats_df$GOM_Plot_Int[1])

GOM_97 <- c("GOM_97", 
            sma_stats_df$GOM_Plot_Elevation[2], 
            sma_stats_df$GOM_Plot_Slope[2], 
            min((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
            max((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
            (sma_stats_df$GOM_Plot_Slope[2]*min((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_df$GOM_Plot_Elevation[2],
            (sma_stats_df$GOM_Plot_Slope[2]*max((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_df$GOM_Plot_Elevation[2],
            sma_stats_df$GOM_R2[2],
            sma_stats_df$GOM_R2_P[2],
            sma_stats_df$GOM_Slope[2],
            sma_stats_df$GOM_Plot_Int[2])

GOM_2 <- c("GOM_2", 
           sma_stats_df$GOM_Plot_Elevation[3], 
           sma_stats_df$GOM_Plot_Slope[3], 
           min((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
           max((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
           (sma_stats_df$GOM_Plot_Slope[3]*min((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_df$GOM_Plot_Elevation[3],
           (sma_stats_df$GOM_Plot_Slope[3]*max((sma_collated_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_df$GOM_Plot_Elevation[3],
           sma_stats_df$GOM_R2[3],
           sma_stats_df$GOM_R2_P[3],
           sma_stats_df$GOM_Slope[3],
           sma_stats_df$GOM_Plot_Int[3])

# this binds the plotting information for all of the regions mean and 95% CI together 
sma_plot_allinfo <- as.data.frame(rbind(GOM_Mean, GOM_97, GOM_2, NYB_Mean, NYB_97, NYB_2, Greenland_Mean, Greenland_97, Greenland_2, Iceland_Mean, Iceland_97, Iceland_2))
colnames(sma_plot_allinfo) <- colnames_plot

# convert to numeric as always !!! 
sma_plot_allinfo <- convert_to_numeric(sma_plot_allinfo, columns = c("elevation", "slope", "min_x","max_x","min_y","max_y","r2","pval","Slope","Int"))

site_colors <- c("New_York_Bight" = "firebrick3", "Iceland" = "lightskyblue", "Gulf_of_Maine" = "goldenrod2", "Historical" = "blue4", "Greenland" = "forestgreen")
sma_site_colors <- c("New_York_Bight" = "firebrick3", "NYB_97" = rgb(0.80,0.15,0.15, 0.50), "NYB_2" = rgb(0.80,0.15,0.15, 0.50), "NYB_Raw" = rgb(0.80,0.15,0.15, 0.35), "NYB_Mean" = rgb(0.80,0.15,0.15,1),
                     "Gulf_of_Maine" = "goldenrod2", "GOM_97" = rgb(0.93,0.71,0.13, 0.50), "GOM_2" = rgb(0.93,0.71,0.13, 0.50), "GOM_Raw" = rgb(0.93,0.71,0.13, 0.35), "GOM_Mean" = rgb(0.93,0.71,0.13, 1),
                     "Iceland" = "lightskyblue", "Iceland_97" = rgb(0.53, 0.81, 0.98, 0.50), "Iceland_2" = rgb(0.53, 0.81, 0.98, 0.50), "Iceland_Raw" = rgb(0.53, 0.81, 0.98, 0.35), "Iceland_Mean" = rgb(0.53, 0.81, 0.98, 1),
                     "Greenland" = "forestgreen", "Greenland_97" = rgb(0.13, 0.55, 0.13, 0.50), "Greenland_2" = rgb(0.13, 0.55, 0.13, 0.50), "Greenland_Raw" = rgb(0.13, 0.55, 0.13, 0.35), "Greenland_Mean" = rgb(0.13, 0.55, 0.13, 1))

mean_animal_df_US <- mean_animal_df
mean_animal_df_US$Foraging_Region <- gsub("New York Bight", "New_York_Bight", mean_animal_df_US$Foraging_Region)
mean_animal_df_US$Foraging_Region <- gsub("Gulf of Maine", "Gulf_of_Maine", mean_animal_df_US$Foraging_Region)


gg_sma_boot_plot <- ggplot(data = mean_animal_df_US, ggplot2::aes(x = LogMean_TL, y = Log_Mean_BV, color= Foraging_Region, fill= Foraging_Region)) +
  ggplot2::geom_point(data = mean_animal_df_US, aes(x = LogMean_TL, y = LogMean_BV, color= Foraging_Region, fill= Foraging_Region), size=2, alpha = 0.5, show.legend = F) +
  ggplot2::geom_segment(data = sma_plot_allinfo ,ggplot2::aes(x= min_x, xend= max_x, y=min_y, yend=max_y, colour= group), inherit.aes = FALSE, size= 1) +
  scale_color_manual(values = sma_site_colors, labels = c("Greenland", "Gulf of Maine", "New York Bight", "Iceland"), breaks = c("Greenland", "Gulf_of_Maine", "New_York_Bight", "Iceland")) +
  scale_alpha(guide = "none") +
  labs(x = "Log (Total Length) (m)", 
       y = bquote('Log (Body Volume) '(m^3))) + 
  theme_bw() +
  coord_cartesian(xlim = c(0.79, 1.22), ylim = c(0.4, 1.83))

gg_sma_boot_plot + theme(legend.position = "bottom", legend.text=element_text(size=19), 
                         legend.title = element_text(size=19),
                         #legend.title=element_text(size=11, face = "bold"),
                         plot.subtitle=element_text(size=22),
                         axis.text = element_text(size=19),
                         axis.title = element_text(size=22),
                         plot.title=element_text(face = "bold", size = 22)) + guides(color = guide_legend(nrow = 2, title = "Foraging Region"))


# Note that shading was editing in using image editing software for clarity, and because geom_ribbon did not like my line segments :( 

# Plotting The Mass SMA!  -------------------------------------------------------------

# SMA Plot 

# This is the same code as the previous section except using the mass data 

colnames_plot <- c("group", "elevation", "slope", "min_x", "max_x", "min_y","max_y", "r2", "pval", "Slope", "Int")

Greenland_Mean <- c("Greenland_Mean", 
                    sma_stats_mass_df$Greenland_Plot_Elevation[1], 
                    sma_stats_mass_df$Greenland_Plot_Slope[1], 
                    min((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                    max((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                    (sma_stats_mass_df$Greenland_Plot_Slope[1]*min((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_mass_df$Greenland_Plot_Elevation[1],
                    (sma_stats_mass_df$Greenland_Plot_Slope[1]*max((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_mass_df$Greenland_Plot_Elevation[1],
                    sma_stats_mass_df$Greenland_R2[1],
                    sma_stats_mass_df$Greenland_R2_P[1],
                    sma_stats_mass_df$Greenland_Slope[1],
                    sma_stats_mass_df$Greenland_Plot_Int[1])

Greenland_97 <- c("Greenland_97", 
                  sma_stats_mass_df$Greenland_Plot_Elevation[2], 
                  sma_stats_mass_df$Greenland_Plot_Slope[2], 
                  min((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                  max((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                  (sma_stats_mass_df$Greenland_Plot_Slope[2]*min((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_mass_df$Greenland_Plot_Elevation[2],
                  (sma_stats_mass_df$Greenland_Plot_Slope[2]*max((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_mass_df$Greenland_Plot_Elevation[2],
                  sma_stats_mass_df$Greenland_R2[2],
                  sma_stats_mass_df$Greenland_R2_P[2],
                  sma_stats_mass_df$Greenland_Slope[2],
                  sma_stats_mass_df$Greenland_Plot_Int[2])

Greenland_2 <- c("Greenland_2", 
                 sma_stats_mass_df$Greenland_Plot_Elevation[3], 
                 sma_stats_mass_df$Greenland_Plot_Slope[3], 
                 min((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                 max((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length),
                 (sma_stats_mass_df$Greenland_Plot_Slope[3]*min((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_mass_df$Greenland_Plot_Elevation[3],
                 (sma_stats_mass_df$Greenland_Plot_Slope[3]*max((mean_mass_data %>% filter(Foraging_Region == "Greenland"))$Total_Length))+sma_stats_mass_df$Greenland_Plot_Elevation[3],
                 sma_stats_mass_df$Greenland_R2[3],
                 sma_stats_mass_df$Greenland_R2_P[3],
                 sma_stats_mass_df$Greenland_Slope[3],
                 sma_stats_mass_df$Greenland_Plot_Int[3])

Iceland_Mean <- c("Iceland_Mean", 
                  sma_stats_mass_df$Iceland_Plot_Elevation[1], 
                  sma_stats_mass_df$Iceland_Plot_Slope[1], 
                  min((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                  max((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                  (sma_stats_mass_df$Iceland_Plot_Slope[1]*min((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_mass_df$Iceland_Plot_Elevation[1],
                  (sma_stats_mass_df$Iceland_Plot_Slope[1]*max((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_mass_df$Iceland_Plot_Elevation[1],
                  sma_stats_mass_df$Iceland_R2[1],
                  sma_stats_mass_df$Iceland_R2_P[1],
                  sma_stats_mass_df$Iceland_Slope[1],
                  sma_stats_mass_df$Iceland_Plot_Int[1])

Iceland_97 <- c("Iceland_97", 
                sma_stats_mass_df$Iceland_Plot_Elevation[2], 
                sma_stats_mass_df$Iceland_Plot_Slope[2], 
                min((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                max((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
                (sma_stats_mass_df$Iceland_Plot_Slope[2]*min((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_mass_df$Iceland_Plot_Elevation[2],
                (sma_stats_mass_df$Iceland_Plot_Slope[2]*max((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_mass_df$Iceland_Plot_Elevation[2],
                sma_stats_mass_df$Iceland_R2[2],
                sma_stats_mass_df$Iceland_R2_P[2],
                sma_stats_mass_df$Iceland_Slope[2],
                sma_stats_mass_df$Iceland_Plot_Int[2])

Iceland_2 <- c("Iceland_2", 
               sma_stats_mass_df$Iceland_Plot_Elevation[3], 
               sma_stats_mass_df$Iceland_Plot_Slope[3], 
               min((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
               max((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length),
               (sma_stats_mass_df$Iceland_Plot_Slope[3]*min((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_mass_df$Iceland_Plot_Elevation[3],
               (sma_stats_mass_df$Iceland_Plot_Slope[3]*max((mean_mass_data %>% filter(Foraging_Region == "Iceland"))$Total_Length))+sma_stats_mass_df$Iceland_Plot_Elevation[3],
               sma_stats_mass_df$Iceland_R2[3],
               sma_stats_mass_df$Iceland_R2_P[3],
               sma_stats_mass_df$Iceland_Slope[3],
               sma_stats_mass_df$Iceland_Plot_Int[3])

NYB_Mean <- c("NYB_Mean", 
              sma_stats_mass_df$NYB_Plot_Elevation[1], 
              sma_stats_mass_df$NYB_Plot_Slope[1], 
              min((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
              max((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
              (sma_stats_mass_df$NYB_Plot_Slope[1]*min((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_mass_df$NYB_Plot_Elevation[1],
              (sma_stats_mass_df$NYB_Plot_Slope[1]*max((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_mass_df$NYB_Plot_Elevation[1],
              sma_stats_mass_df$NYB_R2[1],
              sma_stats_mass_df$NYB_R2_P[1],
              sma_stats_mass_df$NYB_Slope[1],
              sma_stats_mass_df$NYB_Plot_Int[1])

NYB_97 <- c("NYB_97", 
            sma_stats_mass_df$NYB_Plot_Elevation[2], 
            sma_stats_mass_df$NYB_Plot_Slope[2], 
            min((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
            max((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
            (sma_stats_mass_df$NYB_Plot_Slope[2]*min((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_mass_df$NYB_Plot_Elevation[2],
            (sma_stats_mass_df$NYB_Plot_Slope[2]*max((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_mass_df$NYB_Plot_Elevation[2],
            sma_stats_mass_df$NYB_R2[2],
            sma_stats_mass_df$NYB_R2_P[2],
            sma_stats_mass_df$NYB_Slope[2],
            sma_stats_mass_df$NYB_Plot_Int[2])

NYB_2 <- c("NYB_2", 
           sma_stats_mass_df$NYB_Plot_Elevation[3], 
           sma_stats_mass_df$NYB_Plot_Slope[3], 
           min((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
           max((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length),
           (sma_stats_mass_df$NYB_Plot_Slope[3]*min((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_mass_df$NYB_Plot_Elevation[3],
           (sma_stats_mass_df$NYB_Plot_Slope[3]*max((mean_mass_data %>% filter(Foraging_Region == "New_York_Bight"))$Total_Length))+sma_stats_mass_df$NYB_Plot_Elevation[3],
           sma_stats_mass_df$NYB_R2[3],
           sma_stats_mass_df$NYB_R2_P[3],
           sma_stats_mass_df$NYB_Slope[3],
           sma_stats_mass_df$NYB_Plot_Int[3])

GOM_Mean <- c("GOM_Mean", 
              sma_stats_mass_df$GOM_Plot_Elevation[1], 
              sma_stats_mass_df$GOM_Plot_Slope[1], 
              min((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
              max((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
              (sma_stats_mass_df$GOM_Plot_Slope[1]*min((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_mass_df$GOM_Plot_Elevation[1],
              (sma_stats_mass_df$GOM_Plot_Slope[1]*max((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_mass_df$GOM_Plot_Elevation[1],
              sma_stats_mass_df$GOM_R2[1],
              sma_stats_mass_df$GOM_R2_P[1],
              sma_stats_mass_df$GOM_Slope[1],
              sma_stats_mass_df$GOM_Plot_Int[1])

GOM_97 <- c("GOM_97", 
            sma_stats_mass_df$GOM_Plot_Elevation[2], 
            sma_stats_mass_df$GOM_Plot_Slope[2], 
            min((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
            max((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
            (sma_stats_mass_df$GOM_Plot_Slope[2]*min((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_mass_df$GOM_Plot_Elevation[2],
            (sma_stats_mass_df$GOM_Plot_Slope[2]*max((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_mass_df$GOM_Plot_Elevation[2],
            sma_stats_mass_df$GOM_R2[2],
            sma_stats_mass_df$GOM_R2_P[2],
            sma_stats_mass_df$GOM_Slope[2],
            sma_stats_mass_df$GOM_Plot_Int[2])

GOM_2 <- c("GOM_2", 
           sma_stats_mass_df$GOM_Plot_Elevation[3], 
           sma_stats_mass_df$GOM_Plot_Slope[3], 
           min((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
           max((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length),
           (sma_stats_mass_df$GOM_Plot_Slope[3]*min((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_mass_df$GOM_Plot_Elevation[3],
           (sma_stats_mass_df$GOM_Plot_Slope[3]*max((mean_mass_data %>% filter(Foraging_Region == "Gulf_of_Maine"))$Total_Length))+sma_stats_mass_df$GOM_Plot_Elevation[3],
           sma_stats_mass_df$GOM_R2[3],
           sma_stats_mass_df$GOM_R2_P[3],
           sma_stats_mass_df$GOM_Slope[3],
           sma_stats_mass_df$GOM_Plot_Int[3])

sma_plot_allinfo_mass <- as.data.frame(rbind(GOM_Mean, GOM_97, GOM_2, NYB_Mean, NYB_97, NYB_2, Greenland_Mean, Greenland_97, Greenland_2, Iceland_Mean, Iceland_97, Iceland_2))
colnames(sma_plot_allinfo_mass) <- colnames_plot

sma_plot_allinfo_mass <- convert_to_numeric(sma_plot_allinfo_mass, columns = c("elevation", "slope", "min_x","max_x","min_y","max_y","r2","pval","Slope","Int"))

site_colors <- c("New_York_Bight" = "firebrick3", "Iceland" = "lightskyblue", "Gulf_of_Maine" = "goldenrod2", "Historical" = "blue4", "Greenland" = "forestgreen")
sma_site_colors <- c("New_York_Bight" = "firebrick3", "NYB_97" = rgb(0.80,0.15,0.15, 0.50), "NYB_2" = rgb(0.80,0.15,0.15, 0.50), "NYB_Raw" = rgb(0.80,0.15,0.15, 0.35), "NYB_Mean" = rgb(0.80,0.15,0.15,1),
                     "Gulf_of_Maine" = "goldenrod2", "GOM_97" = rgb(0.93,0.71,0.13, 0.50), "GOM_2" = rgb(0.93,0.71,0.13, 0.50), "GOM_Raw" = rgb(0.93,0.71,0.13, 0.35), "GOM_Mean" = rgb(0.93,0.71,0.13, 1),
                     "Iceland" = "lightskyblue", "Iceland_97" = rgb(0.53, 0.81, 0.98, 0.50), "Iceland_2" = rgb(0.53, 0.81, 0.98, 0.50), "Iceland_Raw" = rgb(0.53, 0.81, 0.98, 0.35), "Iceland_Mean" = rgb(0.53, 0.81, 0.98, 1),
                     "Greenland" = "forestgreen", "Greenland_97" = rgb(0.13, 0.55, 0.13, 0.50), "Greenland_2" = rgb(0.13, 0.55, 0.13, 0.50), "Greenland_Raw" = rgb(0.13, 0.55, 0.13, 0.35), "Greenland_Mean" = rgb(0.13, 0.55, 0.13, 1))

mean_animal_df_US <- mean_animal_df
mean_animal_df_US$Foraging_Region <- gsub("New York Bight", "New_York_Bight", mean_animal_df_US$Foraging_Region)
mean_animal_df_US$Foraging_Region <- gsub("Gulf of Maine", "Gulf_of_Maine", mean_animal_df_US$Foraging_Region)
mean_animal_df_US$LogMean_Mass <- log10(mean_animal_df_US$mass)


gg_sma_boot_plot_mass <- ggplot(data = mean_mass_data, ggplot2::aes(x = LogMean_TL, y = Log_Mean_Mass, color= Foraging_Region, fill= Foraging_Region)) +
  ggplot2::geom_point(data = mean_mass_data, aes(x = Total_Length, y = mass, color= Foraging_Region, fill= Foraging_Region), size=2, alpha = 0.5, show.legend = F) +
  ggplot2::geom_segment(data = sma_plot_allinfo_mass ,ggplot2::aes(x= min_x, xend= max_x, y=min_y, yend=max_y, colour= group), inherit.aes = FALSE, size= 1) +
  scale_color_manual(values = sma_site_colors, labels = c("Greenland", "Gulf of Maine", "New York Bight", "Iceland"), breaks = c("Greenland", "Gulf_of_Maine", "New_York_Bight", "Iceland")) +
  scale_alpha(guide = "none") +
  labs(x = "Log (Total Length) (m)", 
       y = bquote('Log (Body Mass) '(kg))) + 
  theme_bw() +
  coord_cartesian(xlim = c(0.81, 1.19), ylim = c(3.5, 4.8))

gg_sma_boot_plot_mass + theme(legend.position = "bottom", legend.text=element_text(size=19), 
                              legend.title = element_text(size=19),
                              #legend.title=element_text(size=11, face = "bold"),
                              plot.subtitle=element_text(size=22),
                              axis.text = element_text(size=19),
                              axis.title = element_text(size=22),
                              plot.title=element_text(face = "bold", size = 22)) + guides(color = guide_legend(nrow = 2, title = "Foraging Region"))


# Note that shading was editing in using image editing software for clarity, and because geom_ribbon did not like my line segments :( 





# Age Class Number Counts -------------------------------

# This section is to calculate the number of adult and juveniles for each region 

age_class_number_func <- function(df){
  temp <- df
  
  nyb_adults <- nrow(temp %>% filter(Foraging_Region == "New York Bight") %>% filter(rep_class == "adult"))
  nyb_juvies <- nrow(temp %>% filter(Foraging_Region == "New York Bight") %>% filter(rep_class == "juvenile"))
  gom_adults <- nrow(temp %>% filter(Foraging_Region == "Gulf of Maine") %>% filter(rep_class == "adult"))
  gom_juvies <- nrow(temp %>% filter(Foraging_Region == "Gulf of Maine") %>% filter(rep_class == "juvenile"))
  ice_adults <- nrow(temp %>% filter(Foraging_Region == "Iceland") %>% filter(rep_class == "adult"))
  ice_juvies <- nrow(temp %>% filter(Foraging_Region == "Iceland") %>% filter(rep_class == "juvenile"))
  gl_adults <- nrow(temp %>% filter(Foraging_Region == "Greenland") %>% filter(rep_class == "adult"))
  gl_juvies <- nrow(temp %>% filter(Foraging_Region == "Greenland") %>% filter(rep_class == "juvenile"))
  
  all_juvies <- nrow(temp %>% filter(rep_class == "juvenile"))
  all_adults <- nrow(temp %>% filter(rep_class == "adult"))
  
  age_results_t <- as.data.frame(rbind(nyb_adults, nyb_juvies, gom_adults, gom_juvies, ice_adults, ice_juvies, gl_adults, gl_juvies, all_juvies, all_adults))
  age_results <- data.table::transpose(age_results_t)
  colnames(age_results) <- rownames(age_results_t)
  
  return(age_results)
  
  
}

age_class <- lapply(pop_boot_all, age_class_number_func)

age_class_mega <- bind_rows(age_class)

# This is a convenient table that calculates age counts for all of the regions 
age_counts <- age_class_mega %>% purrr::map(table)


# This is a plot of the age distribution of all foraging regions 
length_violin <- ggplot(megapop_BCI_df, aes(x = Foraging_Region, y = Total_Length, fill = Foraging_Region, alpha = 0.75, color = Foraging_Region)) + # all regional data, with foraging region as x groups and BCI on y axis. Colors determined by foraging group
  geom_violin(fill = "transparent", draw_quantiles = 0.5, color = "black") + # formatting for violin plot
  geom_violin() + # telling ggplot to make a violin plot
  geom_dotplot(data = mean_animal_df, aes(x = Foraging_Region, y = Total_Length, fill = Foraging_Region, color = Foraging_Region), binaxis = 'y', stackratio = 1, stackgroups = TRUE, stackdir = 'centerwhole', dotsize = 0.12, alpha = 2, binpositions = "bygroup", stroke = NA, method = "histodot") +
  coord_flip() +
  geom_hline(yintercept = 11.47, linetype = "solid", color = "black", size = 0.75) +
  labs(x = "Foraging Region", y = "Total Length") + # labels
  scale_fill_manual(values = site_colors) + # legend 
  scale_color_manual(values = site_colors) +  # legend 
  scale_alpha(guide = "none") +  # legend
  guides(fill="none") +
  guides(color = "none") +
  theme(legend.position = "none") + # title formatting
  theme_bw() + # making theme black and white, removing the default grey background 
  theme(axis.text=element_text(size=12), title=element_text(size=12))
length_violin # plot the plot 


# Mean BCI by Region And Month --------------------------------------------

# This section calculates mean BCI by Region and Month 
megapop_BCI_month <- megapop_BCI_df %>% mutate(month = stringr::str_extract(stringr::str_extract(Date, "\\-\\d{2}\\-"), "\\d{2}"))

# Each region 
megapop_nyb <- megapop_BCI_month %>% filter(Foraging_Region == "New York Bight")
megapop_gom <- megapop_BCI_month %>% filter(Foraging_Region == "Gulf of Maine")
megapop_ice <- megapop_BCI_month %>% filter(Foraging_Region == "Iceland")
megapop_gl <- megapop_BCI_month %>% filter(Foraging_Region == "Greenland")

# BCI By Month
mean_nyb_bci <- mean(megapop_BCI_month %>% filter(Foraging_Region == "New York Bight") %>% pull(BCI))
mean_gl_bci <- mean(megapop_BCI_month %>% filter(Foraging_Region == "Greenland") %>% pull(BCI))
mean_gom_bci <- mean(megapop_BCI_month %>% filter(Foraging_Region == "Gulf of Maine") %>% pull(BCI))
mean_ice_bci <- mean(megapop_BCI_month %>% filter(Foraging_Region == "Iceland") %>% pull(BCI))

# New York Bight BCI Means For Each Month
nyb_bci_means <- list()
for (i in 1:length(unique(megapop_nyb$month))){
  month_uni <- unique(megapop_nyb$month)
  month_uni_i <- month_uni[i]
  nyb_month <- megapop_nyb %>% filter(month == month_uni_i)
  mean_nyb_month <- mean(nyb_month %>% pull(BCI))
  temp <- as.data.frame(rbind(month_uni_i, mean_nyb_month))
  nyb_bci_means[i] <- temp
}

# Gulf of Maine BCI Means For Each Month

gom_bci_means <- list()
for (i in 1:length(unique(megapop_gom$month))){
  month_uni <- unique(megapop_gom$month)
  month_uni_i <- month_uni[i]
  gom_month <- megapop_gom %>% filter(month == month_uni_i)
  mean_gom_month <- mean(gom_month %>% pull(BCI))
  temp <- as.data.frame(rbind(month_uni_i, mean_gom_month))
  gom_bci_means[i] <- temp
}

# Iceland BCI Means For Each Month

ice_bci_means <- list()
for (i in 1:length(unique(megapop_ice$month))){
  month_uni <- unique(megapop_ice$month)
  month_uni_i <- month_uni[i]
  ice_month <- megapop_ice %>% filter(month == month_uni_i)
  mean_ice_month <- mean(ice_month %>% pull(BCI))
  temp <- as.data.frame(rbind(month_uni_i, mean_ice_month))
  ice_bci_means[i] <- temp
}

# Greenland BCI Means For Each Month

gl_bci_means <- list()
for (i in 1:length(unique(megapop_gl$month))){
  month_uni <- unique(megapop_gl$month)
  month_uni_i <- month_uni[i]
  gl_month <- megapop_gl %>% filter(month == month_uni_i)
  mean_gl_month <- mean(gl_month %>% pull(BCI))
  temp <- as.data.frame(rbind(month_uni_i, mean_gl_month))
  gl_bci_means[i] <- temp
}


# SMA Back-calculations  --------------------------------------------------

# This section is how differences in percentage and kgs is calculated from the SMA results 

sma_slope_mean <- sma_stats_mass_df$Greenland_Slope[1] # this is the mean slope for all of the SMA; Greenland is the first / easiest one to pull but they're all the same
sma_slope_low <- sma_stats_mass_df$Greenland_Slope[3] # this is the mean slope for all of the SMA; Greenland is the first / easiest one to pull but they're all the same
sma_slope_high <- sma_stats_mass_df$Greenland_Slope[2] # this is the mean slope for all of the SMA; Greenland is the first / easiest one to pull but they're all the same

sma_gl_int_mean <- sma_stats_mass_df$Greenland_Plot_Int[1] # mean greenland intercept
sma_gom_int_mean <- sma_stats_mass_df$GOM_Plot_Int[1] # mean gom intercept
sma_nyb_int_mean <- sma_stats_mass_df$NYB_Plot_Int[1] # mean nyb intercept
sma_ice_int_mean <- sma_stats_mass_df$Iceland_Plot_Int[1] # mean iceland intercept

sma_gl_int_low <- sma_stats_mass_df$Greenland_Plot_Int[3] # mean greenland intercept
sma_gom_int_low <- sma_stats_mass_df$GOM_Plot_Int[3] # mean gom intercept
sma_nyb_int_low <- sma_stats_mass_df$NYB_Plot_Int[3] # mean nyb intercept
sma_ice_int_low <- sma_stats_mass_df$Iceland_Plot_Int[3] # mean iceland intercept

sma_gl_int_high <- sma_stats_mass_df$Greenland_Plot_Int[2] # mean greenland intercept
sma_gom_int_high <- sma_stats_mass_df$GOM_Plot_Int[2] # mean gom intercept
sma_nyb_int_high <- sma_stats_mass_df$NYB_Plot_Int[2] # mean nyb intercept
sma_ice_int_high <- sma_stats_mass_df$Iceland_Plot_Int[2] # mean iceland intercept


sma_glnybice_mean <- mean(c(sma_gl_int_mean, sma_nyb_int_mean, sma_ice_int_mean)) # combined data for greenland, nyb, and iceland 
sma_glnybice_low <- mean(c(sma_gl_int_low, sma_nyb_int_low, sma_ice_int_low)) # combined data for greenland, nyb, and iceland 
sma_glnybice_high <- mean(c(sma_gl_int_high, sma_nyb_int_high, sma_ice_int_high)) # combined data for greenland, nyb, and iceland 


vertical_shift_mean <- sma_gom_int_mean - sma_glnybice_mean # vertical shift in log log space from mean to mean 
vertical_shift_low <- sma_gom_int_low - sma_glnybice_low # vertical shift in log log space from lower bound to lower bound 
vertical_shift_high <- sma_gom_int_high - sma_glnybice_high # vertical shift in log log space from upper bound to upper bound 


x_values <- c(10, 12) # sample values of whale length 
log_x_values <- log10(x_values) # log 10 values of sample whale lengths 

body_volume_three_mean <- 10^(sma_slope_mean * log_x_values + sma_glnybice_mean) # taking out of log space 
body_volume_gom_mean <- 10^(sma_slope_mean * log_x_values + sma_gom_int_mean) # taking out of log space 

body_volume_three_low <- 10^(sma_slope_low * log_x_values + sma_glnybice_low) # taking out of log space 
body_volume_gom_low <- 10^(sma_slope_low * log_x_values + sma_gom_int_low) # taking out of log space 

body_volume_three_high <- 10^(sma_slope_high * log_x_values + sma_glnybice_high) # taking out of log space 
body_volume_gom_high <- 10^(sma_slope_high * log_x_values + sma_gom_int_high) # taking out of log space 


diff_mean <- body_volume_gom_mean - body_volume_three_mean # difference in mean values 
diff_low <- body_volume_gom_low - body_volume_three_low # difference in lower values 
diff_high <- body_volume_gom_high - body_volume_three_high # difference in higher values 


# Last little plot before we go! 

# Plotting DOY vs BCI trends 

DOYplot_region <- ggplot(mean_animal_df, aes(x = DOY, y = BCI)) +
  geom_point(aes(color = Foraging_Region)) +  # Points colored by Foraging_Region
  scale_color_manual(values = site_colors) +  # Manual scale for colors
  geom_smooth(method = "lm", se = TRUE, color = "black") +     # Regression line for entire dataset
  theme_bw()

# Filter data for only Iceland group
iceland_data <- mean_animal_df[mean_animal_df$Foraging_Region == "Iceland", ]

# Add line for Iceland data
DOYplot_region <- DOYplot_region +
  geom_smooth(data = iceland_data, aes(group = Foraging_Region), method = "lm", se = TRUE, color = "lightskyblue", fill = "lightskyblue")

DOYplot_region




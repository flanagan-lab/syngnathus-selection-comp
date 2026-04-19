############################################################################################
###################### Partitioning the opportunity for selection - ########################
##################### not including unmated individuals after episode 1 ####################
############################################################################################

## objects fem_succ and mal_succ generated for each species in their respective 
## .Rmd documents 

##################################### FEMALES ##############################################
#Create a dataframe to store all of the intermediate values of fitness in
fem_succ_fitness <- data.frame(matrix(ncol = ncol(fem_succ) + 5,
                                      nrow = 0))
colnames(fem_succ_fitness) <- c(colnames(fem_succ),
                                "w1", "W1W2", "w2", "W1W2W3", "w3")

#Create a dataframe to store the final calculations of I and COI in
opp_selection_episodes_fem <- data.frame(matrix(ncol = 16,
                                                nrow = 0))
colnames(opp_selection_episodes_fem) <- c("trial_num", "I_1", "I_2", "coi1_2",
                                          "coi1_2given1", "coi12_2", "coi12_2given1",
                                          "diff_12", "I_12", "I_3", "coi12_3",
                                          "coi12_3given2", "coi123_3", "coi123_3given2",
                                          "diff_123", "I")

#Loop through the individual trials and partition I into all components
for (trial in unique(fem_succ$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- fem_succ[fem_succ$trial_num == trial, ]
  
  #Calculate the absolute pre-copulatory fitness (Eq. 14 Arnold & Wade 1984)
  #This is the same as the calculation of I_s
  tmp$w1 <- tmp$MatingSuccess/mean(tmp$MatingSuccess) #Relative mating success
  p <- (1)/(nrow(tmp) - 1)
  
  I_1 <- sum(p*((tmp$w1 - 1)^2))
  
  #Post-copulatory selection event 1 (Number of eggs transferred per mate) (Eq. 15 Arnold & Wade 1984)
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean).
  tmp$W1W2 <- ifelse(tmp$MatingSuccess == 0,
                     0,
                     tmp$totalEggs/mean(tmp$totalEggs[tmp$MatingSuccess > 0])) #Relative num. eggs transferred
  tmp$w2 <- tmp$W1W2/tmp$w1 #Relative num. eggs transferred PER MATE
  tmp$w2[is.na(tmp$w2)] <- 0
  p1 <- tmp$w1/(nrow(tmp) - 1)
  
  I_2 <- sum(p1*((tmp$w2 - 1)^2))
  
  #Post-copulatory selection event 2 (Proportion eggs developed) (Eq. 16 Arnold & Wade 1984)
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean).
  tmp$W1W2W3 <- ifelse(tmp$MatingSuccess == 0,
                       0,
                       tmp$NumDeveloped/mean(tmp$NumDeveloped[tmp$MatingSuccess > 0])) #Relative num. eggs developed
  tmp$w3 <- tmp$W1W2W3/tmp$W1W2 #Relative PROPORTION eggs developed
  tmp$w3[is.na(tmp$w3)] <- 0
  p2 <- (tmp$w1*tmp$w2)/(nrow(tmp)-1)
  
  I_3 <- sum(p2*((tmp$w3 - 1)^2))
  
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean and variance).
  I_12 <- var(tmp$W1W2[tmp$W1W2 > 0])/(mean(tmp$W1W2[tmp$W1W2 > 0])^2)
  
  #Total opportunity for selection
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean and variance).
  I <- var(tmp$W1W2W3[tmp$W1W2W3 > 0])/(mean(tmp$W1W2W3[tmp$W1W2W3 > 0])^2)
  
  ##Calculating the co-variances
  coi1_2 <- cov(tmp$w1,tmp$w2) ##COI(1, 2) Eq. A1
  coi1_2given1 <- ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*tmp$w2)) - ((1/(nrow(tmp) - 1))*sum(tmp$w1^2)) ##COI(1, 2|1) Eq. A2
  
  coi12_2 <- (1/(nrow(tmp) - 1))*sum(tmp$w1*(tmp$w2^2)) - (1/(nrow(tmp) - 1))*sum(tmp$w2) ##COI(12, 2) Eq. A3
  coi12_2given1 <- ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*(tmp$w2^2))) - ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*tmp$w2)) ##COI(12, 2|1) Eq. A4
  
  diff_12 <- coi12_2given1 - coi12_2 ##Used in Eq. A5 for calculating I_12
  
  coi12_3 <- sum(p*((tmp$w1*tmp$w2) - 1)*(tmp$w3 - 1)) ##COI(12, 3) Eq. A6
  coi12_3given2 <- ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*(tmp$w2^2)*(tmp$w3))) - ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*(tmp$w2^2))) ##COI(12, 3|2) Eq. A9
  
  coi123_3 <- sum(p*((tmp$w1*tmp$w2*tmp$w3) - 1)*(tmp$w3 - 1)) ##COI(123, 3) Eq. A10
  coi123_3given2 <- sum(p2*(((tmp$w1*tmp$w2*tmp$w3) - 1)*(tmp$w3 - 1))) ##COI(123, 3|2) Eq. All
  
  diff_123 <- coi123_3given2 - coi123_3 ##Used in Eq. A12 for calculating I
  
  
  #Combining all of the selection values (Is) and covariances (COIs) and saving the output
  ##Use for S. scovelli dataset:
  #trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
  #                             trial))
  trial_num <- trial
  selection <- cbind(trial_num, I_1, I_2, coi1_2, coi1_2given1, coi12_2, coi12_2given1,
                     diff_12, I_12, I_3, coi12_3, coi12_3given2, coi123_3, coi123_3given2,
                     diff_123, I)
  
  opp_selection_episodes_fem <- rbind(opp_selection_episodes_fem, selection)
  
  #Save the intermediate values
  fem_succ_fitness <- rbind(fem_succ_fitness, tmp)
}

####################################### MALES ##############################################
#Create a dataframe to store all of the intermediate values of fitness in
mal_succ_fitness <- data.frame(matrix(ncol = ncol(mal_succSS) + 5,
                                      nrow = 0))
colnames(mal_succ_fitness) <- c(colnames(mal_succSS),
                                "w1", "W1W2", "w2", "W1W2W3", "w3")

#Create a dataframe to store the final calculations of I in
opp_selection_episodes_mal <- data.frame(matrix(ncol = 16,
                                                nrow = 0))
colnames(opp_selection_episodes_mal) <- c("trial_num", "I_1", "I_2", "coi1_2",
                                          "coi1_2given1", "coi12_2", "coi12_2given1",
                                          "diff_12", "I_12", "I_3", "coi12_3",
                                          "coi12_3given2", "coi123_3", "coi123_3given2",
                                          "diff_123", "I")

for (trial in unique(mal_succ$trial_num)) {
  
  #Subset the overall dataframe to work with an individual trial
  tmp <- mal_succ[mal_succ$trial_num == trial, ]
  
  #Calculate the absolute pre-copultory fitness (Eq. 14 Arnold & Wade 1984)
  tmp$w1 <- tmp$MatingSuccess/mean(tmp$MatingSuccess) #Relative mating success
  p <- (1)/(nrow(tmp) - 1)
  
  I_1 <- sum(p*((tmp$w1 - 1)^2))
  
  #Post-copulatory selection event 1 (Number of eggs transferred per mate) (Eq. 15 Arnold & Wade 1984)
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean).
  tmp$W1W2 <- ifelse(tmp$MatingSuccess == 0,
                     0,
                     tmp$totalEggs/mean(tmp$totalEggs[tmp$MatingSuccess > 0])) #Relative num. eggs transferred
  tmp$w2 <- tmp$W1W2/tmp$w1 #Relative num. eggs transferred PER MATE
  tmp$w2[is.na(tmp$w2)] <- 0
  p1 <- tmp$w1/(nrow(tmp) - 1)
  
  I_2 <- sum(p1*((tmp$w2 - 1)^2))
  
  #Post-copulatory selection event 2 (Proportion eggs developed) (Eq. 16 Arnold & Wade 1984)
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean).
  tmp$W1W2W3 <- ifelse(tmp$MatingSuccess == 0,
                       0,
                       tmp$NumDeveloped/mean(tmp$NumDeveloped[tmp$MatingSuccess > 0])) #Relative num. eggs developed
  tmp$w3 <- tmp$W1W2W3/tmp$W1W2 #Relative PROPORTION eggs developed
  tmp$w3[is.na(tmp$w3)] <- 0
  p2 <- (tmp$w1*tmp$w2)/(nrow(tmp)-1)
  
  I_3 <- sum(p2*((tmp$w3 - 1)^2))
  
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean and variance).
  I_12 <- var(tmp$W1W2[tmp$W1W2 > 0])/(mean(tmp$W1W2[tmp$W1W2 > 0])^2)
  
  #Total opportunity for selection
  ##Unmated individuals are not included in the generation of relative fitness (i.e, removed
  ##from the calculation of the mean and variance).
  I <- var(tmp$W1W2W3[tmp$W1W2W3 > 0])/(mean(tmp$W1W2W3[tmp$W1W2W3 > 0])^2)
  
  ##Calculating the co-variances
  coi1_2 <- cov(tmp$w1,tmp$w2) ##COI(1, 2) Eq. A1
  coi1_2given1 <- ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*tmp$w2)) - ((1/(nrow(tmp) - 1))*sum(tmp$w1^2)) ##COI(1, 2|1) Eq. A2
  
  coi12_2 <- (1/(nrow(tmp) - 1))*sum(tmp$w1*(tmp$w2^2)) - (1/(nrow(tmp) - 1))*sum(tmp$w2) ##COI(12, 2) Eq. A3
  coi12_2given1 <- ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*(tmp$w2^2))) - ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*tmp$w2)) ##COI(12, 2|1) Eq. A4
  
  diff_12 <- coi12_2given1 - coi12_2 ##Used in Eq. A5 for calculating I_12
  
  coi12_3 <- sum(p*((tmp$w1*tmp$w2) - 1)*(tmp$w3 - 1)) ##COI(12, 3) Eq. A6
  coi12_3given2 <- ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*(tmp$w2^2)*(tmp$w3))) - ((1/(nrow(tmp) - 1))*sum((tmp$w1^2)*(tmp$w2^2))) ##COI(12, 3|2) Eq. A9
  
  coi123_3 <- sum(p*((tmp$w1*tmp$w2*tmp$w3) - 1)*(tmp$w3 - 1)) ##COI(123, 3) Eq. A10
  coi123_3given2 <- sum(p2*(((tmp$w1*tmp$w2*tmp$w3) - 1)*(tmp$w3 - 1))) ##COI(123, 3|2) Eq. All
  
  diff_123 <- coi123_3given2 - coi123_3 ##Used in Eq. A12 for calculating I
  
  
  #Combining all of the selection values (Is) and covariances (COIs) and saving the output
  ##Use for S. scovelli dataset:
  #trial_num <- as.numeric(gsub("^(C)(\\d)", "\\2",
  #                             trial))
  trial_num <- trial
  selection <- cbind(trial_num, I_1, I_2, coi1_2, coi1_2given1, coi12_2, coi12_2given1,
                     diff_12, I_12, I_3, coi12_3, coi12_3given2, coi123_3, coi123_3given2,
                     diff_123, I)
  
  opp_selection_episodes_mal <- rbind(opp_selection_episodes_mal, selection)
  
  #Save the intermediate values
  mal_succ_fitness <- rbind(mal_succ_fitness, tmp)
}

####################### SUMMARY STATS #######################################################
#Merge the selection coefficients from males and females into one dataset to 
#make life easier
opp_selection_episodes_fem$Sex <- "F"
opp_selection_episodes_mal$Sex <- "M"

opp_selection_episodes_all <- rbind(opp_selection_episodes_fem, opp_selection_episodes_mal)

#List the columns of interest
columns <- c("I_1", "I_2", "coi1_2", "coi1_2given1", "coi12_2", "coi12_2given1", "diff_12",
             "I_12", "I_3", "coi12_3", "coi12_3given2", "coi123_3", "coi123_3given2", 
             "diff_123", "I")

#Create a dataframe to store the final values in
opp_episodes_average <- data.frame(matrix(ncol = 5,
                                          nrow = 0))
colnames(opp_episodes_average) <- c("Average", "Variance", "Interval", 
                                    "Episode_sel", "Sex")

#Calculate the critical value
crit <- qt(p = 0.975, df = (nrow(opp_selection_episodes_fem) - 1))

for (j in 1:length(columns)) {
  
  col_name <- columns[[j]]
  
  #Calculate the means
  mean <- t(t(tapply(opp_selection_episodes_all[, colnames(opp_selection_episodes_all) 
                                                == col_name], 
                     opp_selection_episodes_all$Sex, 
                     mean, na.rm = TRUE)))
  
  #Calculate the variance
  var <- t(t(tapply(opp_selection_episodes_all[, colnames(opp_selection_episodes_all) 
                                               == col_name], 
                    opp_selection_episodes_all$Sex, 
                    var, na.rm = TRUE)))
  
  #Calculate standard error
  se <- t(t(tapply(opp_selection_episodes_all[, colnames(opp_selection_episodes_all) 
                                              == col_name], 
                   opp_selection_episodes_all$Sex, 
                   function(x){
                     sqrt(var(x))/sqrt(length(x))
                   })))
  
  #Calculate the value that is added and subtracted from the mean
  int <- se*crit
  
  #Combine the data together
  episode <- as.data.frame(cbind(mean, var, int))
  colnames(episode) <- c("Average", "Variance", "Interval")
  
  episode$Episode_sel <- col_name
  episode$Sex <- rownames(episode)
  
  rownames(episode) <- NULL
  
  opp_episodes_average <- rbind(opp_episodes_average, episode)
  
}

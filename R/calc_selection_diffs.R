# function to calculate selection differentials
# intended for use on a single population (trial in this case)
# reformatted version of Coley Tosto's code in selection_analysis_floridae.Rmd
calc_selection_diffs<-function(reprod_succ){

  # calculate fitness due to mating success (relative mating)
  reprod_succ$fit1 <- reprod_succ$MatingSuccess/mean(reprod_succ$MatingSuccess)
  
  # calculate eggs per mate
  reprod_succ$eggs_per_mate <- reprod_succ$totalEggs/reprod_succ$MatingSuccess
  ##If mating success = 0, eggs_per_mate = NA and it not included in the calculation
  ##of the relative fitness moving forward
  reprod_succ$fit2 <- ifelse(reprod_succ$MatingSuccess > 0,
                     reprod_succ$eggs_per_mate/mean(reprod_succ$eggs_per_mate, na.rm = TRUE),
                     0) #Relative eggs transferred
  
  #Calculate fitness relating to post-mating selection (eggs that developed)
  reprod_succ$prop_dev <- (reprod_succ$NumDeveloped/reprod_succ$MatingSuccess)/reprod_succ$eggs_per_mate
  reprod_succ$fit3 <- ifelse(reprod_succ$MatingSuccess > 0,
                     reprod_succ$prop_dev/mean(reprod_succ$prop_dev, na.rm = TRUE),
                     0)
  
  #Standardizing the trait value to have a mean of 0 and sd of unity
  reprod_succ$StdLength <- (reprod_succ$svl - mean(reprod_succ$svl))/sd(reprod_succ$svl)
  
  #Calculating the absolute selection differentials (s)
  s1 <- cov(reprod_succ$svl, reprod_succ$fit1)
  s12 <- cov(reprod_succ$svl, reprod_succ$fit2)
  s123 <- cov(reprod_succ$svl, reprod_succ$fit3)
  s2 <- s12 - s1
  s3 <- s123 - s12
  
  #Calculating the standardized selection differentials (s')
  s1_prime <- cov(reprod_succ$StdLength, reprod_succ$fit1)
  s12_prime <- cov(reprod_succ$StdLength, reprod_succ$fit2)
  s123_prime <- cov(reprod_succ$StdLength, reprod_succ$fit3)
  s2_prime <- s12_prime - s1_prime
  s3_prime <- s123_prime - s12_prime
  
  #Combining all of the selection differentials (s, s') and saving the output
  selection <- cbind(s1, s2, s3, s12, s123, 
                     s1_prime, s2_prime, s3_prime, s12_prime, s123_prime)
  return(selection)
}

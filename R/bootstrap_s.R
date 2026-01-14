# function to bootstrap the calculation of s and s'
bootstrap_s <- function(num_bootstraps,
                        fem_succ, 
                        mal_succ){
  # define the partitions
  all_partitions<-c("s1", "s2", "s3", "s12", "s123", 
  "s1_prime", "s2_prime", "s3_prime", "s12_prime", "s123_prime")
  
  boot_means<- data.frame(matrix(ncol = length(all_partitions)*2,
                                 nrow = num_bootstraps,
                                 data = NA))
  colnames(boot_means)<-c(paste0(all_partitions,"_fem"),
                          paste0(all_partitions,"_mal"))
  
  for(i in 1:num_bootstraps){
    # set up a data frames for storing each boostrap
    this_boot_fem<- data.frame(matrix(ncol = ncol(fem_succ),
                                      nrow = length(unique(fem_succ$trial_num))))
    colnames(this_boot_fem)<-colnames(fem_succ)
    rownames(this_boot_fem)<-paste("trial",unique(fem_succ$trial_num))
    
    this_boot_mal<- data.frame(matrix(ncol = ncol(mal_succ),
                                      nrow = length(unique(mal_succ$trial_num))))
    colnames(this_boot_mal)<-colnames(mal_succ)
    rownames(this_boot_mal)<-paste("trial",unique(mal_succ$trial_num))
   
    count <- 1
    # Resample within each trial
    for(t in unique(fem_succ$trial_num)){
      # females
      this_trial<-fem_succ[fem_succ$trial_num %in% t,]
      this_boot_fem[count,] <- this_trial[sample(1:nrow(this_trial), 1, replace=TRUE),]
      
      # males
      this_trial<-mal_succ[mal_succ$trial_num %in% t,]
      this_boot_mal[count,] <- this_trial[sample(1:nrow(this_trial), 1, replace=TRUE),]
      
      count <- count + 1
    }
    
    resamp_sf <- calc_selection_diffs(this_boot_fem)
    resamp_sm <- calc_selection_diffs(this_boot_mal)
    
    # calculate means of partitioned I 
    boot_means[i,]<-c(resamp_sf,
                      resamp_sm)
    if(!complete.cases(boot_means[i,])){
      print(paste("bootstrap rep",i,"had samples with no mating success"))
      boot_means[i,][is.na(boot_means[i,])]<-0
    }
  }
  return(boot_means)
}

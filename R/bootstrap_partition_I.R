# bootstrap function

bootstrap_partition_I <- function(num_bootstraps,
                                  fem_succ,
                                  mal_succ){
  # define the partitions
  all_partitions<-c("I_1", "I_2", "coi1_2",
                    "coi1_2given1", "coi12_2", "coi12_2given1",
                    "diff_12", "I_12", "I_3", "coi12_3",
                    "coi12_3given2", "coi123_3", "coi123_3given2",
                    "diff_123", "I")
  
  boot_means<- data.frame(matrix(ncol = length(all_partitions)*2,
                                 nrow = num_bootstraps,
                                 data = NA))
  colnames(boot_means)<-c(paste0(all_partitions,"_fem"),
                          paste0(all_partitions,"_mal"))
  
  for(i in 1:num_bootstraps){
    
    # set up a data frames for storing each boostrap
    this_boot_fem<- data.frame(matrix(ncol = length(all_partitions),
                                      nrow = length(unique(fem_succ$trial_num))))
    colnames(this_boot_fem)<-all_partitions
    rownames(this_boot_fem)<-paste("trial",unique(fem_succ$trial_num))
    
    this_boot_mal<- data.frame(matrix(ncol = length(all_partitions),
                                      nrow = length(unique(mal_succ$trial_num))))
    colnames(this_boot_mal)<-all_partitions
    rownames(this_boot_mal)<-paste("trial",unique(mal_succ$trial_num))
    
    count <- 1
    # Resample within each trial
    for(t in unique(fem_succ$trial_num)){
      # females
      this_trial<-fem_succ[fem_succ$trial_num %in% t,]
      resamp <- this_trial[sample(1:nrow(this_trial), nrow(this_trial), replace=TRUE),]
      resamp_If <- partition_I(resamp)
      
      this_boot_fem[count,]<- resamp_If
      
      
      # males
      this_trial<-mal_succ[mal_succ$trial_num %in% t,]
      resamp <- this_trial[sample(1:nrow(this_trial), nrow(this_trial), replace=TRUE),]
      resamp_Im <- partition_I(resamp)
      this_boot_mal[count,]<-resamp_Im
      
      count <- count + 1
    }
    
    # calculate means of partitioned I 
    boot_means[i,]<-c(colMeans(this_boot_fem),
                      colMeans(this_boot_mal))
    if(!complete.cases(boot_means[i,])){
      browser()
    }
  }
  return(boot_means)
}
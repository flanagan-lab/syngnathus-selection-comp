# function to partition the opportunity for selection for pipefish
partition_I <- function(reprod_succ){

    #Calculate the absolute pre-copulatory fitness (Eq. 14 Arnold & Wade 1984)
    #This is the same as the calculation of I_s
    reprod_succ$w1 <- reprod_succ$MatingSuccess/mean(reprod_succ$MatingSuccess) #Relative mating success
    p <- (1)/(nrow(reprod_succ) - 1)
    
    I_1 <- sum(p*((reprod_succ$w1 - 1)^2))
    
    #Post-copulatory selection event 1 (Number of eggs transferred per mate) (Eq. 15 Arnold & Wade 1984)
    reprod_succ$W1W2 <- reprod_succ$totalEggs/mean(reprod_succ$totalEggs) #Relative num. eggs transferred
    reprod_succ$w2 <- reprod_succ$W1W2/reprod_succ$w1 #Relative num. eggs transferred PER MATE
    reprod_succ$w2[is.na(reprod_succ$w2)] <- 0
    p1 <- reprod_succ$w1/(nrow(reprod_succ) - 1)
    
    I_2 <- sum(p1*((reprod_succ$w2 - 1)^2))
    
    #Post-copulatory selection event 2 (Proportion eggs developed) (Eq. 16 Arnold & Wade 1984)
    reprod_succ$W1W2W3 <- reprod_succ$NumDeveloped/mean(reprod_succ$NumDeveloped) #Relative num. eggs developed
    reprod_succ$w3 <- reprod_succ$W1W2W3/reprod_succ$W1W2 #Relative PROPORTION eggs developed
    reprod_succ$w3[is.na(reprod_succ$w3)] <- 0
    p2 <- (reprod_succ$w1*reprod_succ$w2)/(nrow(reprod_succ)-1)
    
    I_3 <- sum(p2*((reprod_succ$w3 - 1)^2))
    
    I_12 <- var(reprod_succ$W1W2)/(mean(reprod_succ$W1W2)^2)
    
    #Total opportunity for selection
    I <- var(reprod_succ$W1W2W3)/(mean(reprod_succ$W1W2W3)^2)
    
    ##Calculating the co-variances
    coi1_2 <- cov(reprod_succ$w1,reprod_succ$w2) ##COI(1, 2) Eq. A1
    coi1_2given1 <- ((1/(nrow(reprod_succ) - 1))*sum((reprod_succ$w1^2)*reprod_succ$w2)) - ((1/(nrow(reprod_succ) - 1))*sum(reprod_succ$w1^2)) ##COI(1, 2|1) Eq. A2
    
    coi12_2 <- (1/(nrow(reprod_succ) - 1))*sum(reprod_succ$w1*(reprod_succ$w2^2)) - (1/(nrow(reprod_succ) - 1))*sum(reprod_succ$w2) ##COI(12, 2) Eq. A3
    coi12_2given1 <- ((1/(nrow(reprod_succ) - 1))*sum((reprod_succ$w1^2)*(reprod_succ$w2^2))) - ((1/(nrow(reprod_succ) - 1))*sum((reprod_succ$w1^2)*reprod_succ$w2)) ##COI(12, 2|1) Eq. A4
    
    diff_12 <- coi12_2given1 - coi12_2 ##Used in Eq. A5 for calculating I_12
    
    coi12_3 <- sum(p*((reprod_succ$w1*reprod_succ$w2) - 1)*(reprod_succ$w3 - 1)) ##COI(12, 3) Eq. A6
    coi12_3given2 <- ((1/(nrow(reprod_succ) - 1))*sum((reprod_succ$w1^2)*(reprod_succ$w2^2)*(reprod_succ$w3))) - ((1/(nrow(reprod_succ) - 1))*sum((reprod_succ$w1^2)*(reprod_succ$w2^2))) ##COI(12, 3|2) Eq. A9
    
    coi123_3 <- sum(p*((reprod_succ$w1*reprod_succ$w2*reprod_succ$w3) - 1)*(reprod_succ$w3 - 1)) ##COI(123, 3) Eq. A10
    coi123_3given2 <- sum(p2*(((reprod_succ$w1*reprod_succ$w2*reprod_succ$w3) - 1)*(reprod_succ$w3 - 1))) ##COI(123, 3|2) Eq. All
    
    diff_123 <- coi123_3given2 - coi123_3 ##Used in Eq. A12 for calculating I
    
    
    #Combining all of the selection values (Is) and covariances (COIs) and saving the output
    selection <- cbind(I_1, I_2, coi1_2, coi1_2given1, coi12_2, coi12_2given1,
                       diff_12, I_12, I_3, coi12_3, coi12_3given2, coi123_3, coi123_3given2,
                       diff_123, I)
    colnames(selection)<-  c("I_1", "I_2", "coi1_2",
                             "coi1_2given1", "coi12_2", "coi12_2given1",
                             "diff_12", "I_12", "I_3", "coi12_3",
                             "coi12_3given2", "coi123_3", "coi123_3given2",
                             "diff_123", "I")
    return(selection)
  
}
#Compute power from relevant inputs

ComputePowerBayesianGeneral <- function(sample_size=500,
                                 response_prob = c(0.5,0.9,0.7,0.2,0.2,0.8,0.2,0.7),
                                 stage_one_trt_one_response_prob = 0.7,
                                 stage_one_trt_two_response_prob = 0.5,
                                 rejection_indices=c(1,2,4,5,6)){
  
  #Arguments
  #sample_size: total sample size in SMART study
  #response_prob: probability of response for each of eight embedded treatment sequences
  #stage_one_trt_one_response_prob: probability of response at the end of stage-1 to treatment a1=1
  #stage_one_trt_two_response_prob: probability of response at the end of stage-1 to treatment a1=0
  #rejection_indices: indices of embedded DTRs to exclude from the set of best in power calculation

  
  upper_limit <- array(NA, dim=c(1000,10,8))



  #Stage-1 response indicator (1 if responder, 0 otherwise)
  s <- rep(NA,sample_size)
  
  
  for (i in 1:1000) {
  
  #Stage-1 treatment indicator
  a1 <- rbinom(sample_size,1,0.5)
  s[a1==1] <- rbinom(length(which(a1==1)),1,stage_one_trt_two_response_prob)
  s[a1==0] <- rbinom(length(which(a1==0)),1,stage_one_trt_one_response_prob)

  #Stage-2 randomization indicator for responders to stage-1 treatment
  a2r<-2*rbinom(sample_size,1,0.5)-1
  
  #Stage-2 randomization indicator for non-responders to stage-1 treatment
  a2nr <- 2*rbinom(sample_size,1,0.5)-1
  
  
  #End of study response indicators for each of the eight embedded treatment sequences (not DTRs) (1 is response, 0 is non-response)
  y1 <- rbinom(length(which(a1==1&s==1&a2r==1)),1,mean(response_prob[1]))
  y2 <- rbinom(length(which(a1==1&s==1&a2r==-1)),1,mean(response_prob[2]))
  y3 <- rbinom(length(which(a1==1&s==0&a2nr==1)),1,mean(response_prob[3]))
  y4 <- rbinom(length(which(a1==1&s==0&a2nr==-1)),1,mean(response_prob[4]))
  y5 <- rbinom(length(which(a1==0&s==1&a2r==1)),1,mean(response_prob[5]))
  y6 <- rbinom(length(which(a1==0&s==1&a2r==-1)),1,mean(response_prob[6]))
  y7 <- rbinom(length(which(a1==0&s==0&a2nr==1)),1,mean(response_prob[7]))
  y8 <- rbinom(length(which(a1==0&s==0&a2nr==-1)),1,mean(response_prob[8]))
  
  
  for (j in 1:10) {
  #Draw 100 draws from the posterior of the probability of response at the end of the study for each of the eight embedded treatment sequences.
  p_1_results <- rbeta(1000, shape1 = sum(y1)+1, shape2=sum(a1==1&s==1&a2r==1)-sum(y1)+1)
  p_2_results <- rbeta(1000, shape1 = sum(y2)+1,shape2 = sum(a1==1&s==1&a2r==-1)-sum(y2)+1)
  p_3_results <- rbeta(1000, shape1 = sum(y3)+1,shape2 = sum(a1==1&s==0&a2nr==1)-sum(y3)+1)
  p_4_results <- rbeta(1000, shape1 = sum(y4)+1,shape2 = sum(a1==1&s==0&a2nr==-1)-sum(y4)+1)
  p_5_results <- rbeta(1000, shape1 = sum(y5)+1,shape2 = sum(a1==0&s==1&a2r==1)-sum(y5)+1)
  p_6_results <- rbeta(1000, shape1 = sum(y6)+1,shape2 = sum(a1==0&s==1&a2r==-1)-sum(y6)+1)
  p_7_results <- rbeta(1000, shape1 = sum(y7)+1, shape2 = sum(a1==0&s==0&a2nr==1)-sum(y7)+1)
  p_8_results <- rbeta(1000, shape1 = sum(y8)+1, shape2 = sum(a1==0&s==0&a2nr==-1)-sum(y8)+1)
  


  #Draw 1000 draws from the posterior of the probability of response at the end of stage 1 for stage-1 treatment 1 and 0, respectively.
  s1 <- rbeta(1000, sum(s[a1==1])+1,sum(a1==1)-sum(s[a1==1])+1)
  s2 <- rbeta(1000, sum(s[a1==0])+1,sum(a1==0)-sum(s[a1==0])+1)
  
  
  #Compute embedded DTR end of study response probability draws using Robin's G-computation formula
  thetadraws <- cbind(p_1_results*(s1)+p_3_results*(1-(s1)),
        p_1_results*(s1)+p_4_results*(1-(s1)),
        p_2_results*(s1)+p_3_results*(1-(s1)),
        p_2_results*(s1)+p_4_results*(1-(s1)),
        p_5_results*(s2)+p_7_results*(1-(s2)),
        p_5_results*(s2)+p_8_results*(1-(s2)),
        p_6_results*(s2)+p_7_results*(1-(s2)),
        p_6_results*(s2)+p_8_results*(1-(s2)))
  
  
  #Perform Bayesian MCB
  thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
  
  log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=1000,ncol=8)
  
  
  rank_matrix <- apply(log_OR_matrix,2,rank,ties.method = 'min')
  
  rank_max <- apply(rank_matrix,1,max)
  
  
  new_dat <- apply(log_OR_matrix[,],2,sort)
  
  ranks_quantile <- floor(quantile(rank_max,1-0.05))
  
  upper_limit[i,j,] <-new_dat[ranks_quantile,]
  print(i)
  }
  
  }
  if (length(rejection_indices)==1) {
    return(mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0))
    
  } else {
    return(mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod)))
  }
}





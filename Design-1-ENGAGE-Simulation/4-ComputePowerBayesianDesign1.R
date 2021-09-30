#This function estimates the power given parameters about a SMART study of design-1/ENGAGE-type SMART.

ComputePowerBayesianDesign1 <- function(sample_size=100,
                                 response_prob=c(0.5,0.9,0.3,0.7,0.5,0.8),
                                 stage_one_trt_one_response_prob = 0.7,
                                 stage_one_trt_two_response_prob = 0.5,
                                 threshold,
                                 alpha = 0.05)
                                 {

  # Arguments:
  # sample_size: total SMART study sample size
  # response_prob: probability of response for each of six embedded treatment sequences
  # stage_one_trt_one_response_prob: probability of response to stage-1 treatment given initial treatment one
  # stage_one_trt_two_response_prob: probability of response to stage-1 treatment given initial treatment two
  # threshold: threshold for exclusion from set of best
  # alpha: probability of excluding best EDTR from set of best
  
  pb <- txtProgressBar(min=0,
                       max=1000,
                       initial = 0,
                       style=3)
  
  upper_limit <- array(NA, dim=c(1000,10,4))

for (i in 1:1000) {
  
  #Stage-1 treatment indicator
  a1_binom <-rbinom(sample_size,1,0.5)
  
  #Stage-2 randomization indicator given stage-1 treatment was 1
  a21_binom <-rbinom(sample_size,1,0.5)
  #Stage-2 randomization indicator given stage-1 treatment was 0
  a22_binom <-rbinom(sample_size,1,0.5)
  
  #Posterior probability of being randomized to 1 vs 0 in stage 1
  a1 <- rbeta(2000, sum(a1_binom)+1,sample_size-sum(a1_binom)+1)
  
  #Posterior probability of being randomized to 1 vs. 0 in stage-2
  #This is for first stage being 1
  a21 <- rbeta(2000, sum(a21_binom)+1,sample_size-sum(a21_binom)+1)
  #This is for first stage being 0
  a22 <- rbeta(2000, sum(a22_binom)+1,sample_size-sum(a22_binom)+1)
  
  #Indicator of response to stage 1 treatment 1
  s1_binom <- rbinom(floor(sample_size*mean(a1_binom)),1,stage_one_trt_one_response_prob)
  #Posterior probability of response to stage-1 treatment 1
  s1 <- rbeta(2000, sum(s1_binom)+1,sum(a1_binom)-sum(s1_binom)+1)
  
  #Indicator of response to stage-1 treatment 0
  s2_binom <- rbinom(floor(sample_size*mean(1-a1_binom)),1,stage_one_trt_two_response_prob)
  
  #Posterior probability of response to stage-1 treatment 0
  s2 <- rbeta(2000, sum(s2_binom)+1,sum((1-a1_binom))-sum(s2_binom)+1)

  #Response indicator at end of study for each of the embedded treatment sequences.
  y_1_results <- rbinom(ceiling(mean((sample_size*s1*a1))), 1,response_prob[1])
  y_2_results <- rbinom(ceiling(mean((sample_size*a21*(1-s1)*a1))),1,response_prob[2])
  y_3_results <- rbinom(ceiling(mean((sample_size*(1-a21)*(1-s1)*a1))),1,response_prob[3])
  
  y_4_results <- rbinom(ceiling(mean((sample_size*s2*(1-a1)))),1,response_prob[4])
  y_5_results <- rbinom(ceiling(mean((sample_size*a22*(1-s2)*(1-a1)))),1,response_prob[5])
  y_6_results <- rbinom(ceiling(mean((sample_size*(1-a22)*(1-s2)*(1-a1)))),1,response_prob[6])
  
  for (j in 1:10) {
    
  #M=1000, number of MC samples
  #j is number of times drawing the phi's
  #i is the number of datasets
  
  #Posterior probability of response for each of the six embedded treatment sequences
  p_1_results <- rbeta(2000, shape1 = sum(y_1_results)+1,length(y_1_results)-sum(y_1_results)+1)
  p_2_results <- rbeta(2000, shape1 = sum(y_2_results)+1,length(y_2_results)-sum(y_2_results)+1)
  p_3_results <- rbeta(2000, shape1 = sum(y_3_results)+1,length(y_3_results)-sum(y_3_results)+1)
  
  p_4_results <- rbeta(2000, shape1 = sum(y_4_results)+1,length(y_4_results)-sum(y_4_results)+1)
  p_5_results <- rbeta(2000, shape1 = sum(y_5_results)+1,length(y_5_results)-sum(y_5_results)+1)
  p_6_results <- rbeta(2000, shape1 = sum(y_6_results)+1,length(y_6_results)-sum(y_6_results)+1)
  
  
  #Transform draws from treatment sequence response probabilities and stage-1 treatment response probabilities to 
  #embedded DTR response probabilities using Robins' G-computation method
  thetadraws <- cbind(p_1_results*(s1)+p_2_results*((1-s1)),
                      p_1_results*(s1)+p_3_results*((1-s1)),
                      p_4_results*(s2)+p_5_results*((1-s2)),
                      p_4_results*(s2)+p_6_results*((1-s2)))
  
  #Compute set of best
  
    logORThreshold <- LogOR(response_prob,
                            stage_one_trt_one_response_prob,
                            stage_one_trt_two_response_prob)
    
  thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
  
  log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=2000,ncol=4)
  
  
  rank_matrix <- apply(log_OR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
  
  rank_max <- apply(rank_matrix,1,max)
  
  
  new_dat <- apply(log_OR_matrix[,],2,sort)
  
  ranks_quantile <- ceiling(quantile(rank_max,1-alpha))
  
  upper_limit[i,j,] <-new_dat[ranks_quantile,]
  
  
  
  }
  

  

  setTxtProgressBar(pb, i)
  
}
  close(pb)
  
    rejection_indices <- which(abs(logORThreshold)>threshold)

  if (length(rejection_indices)==1){
    return(mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0))
    
  } else {
    return(mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod)))
  }

}



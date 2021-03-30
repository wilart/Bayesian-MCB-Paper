##Simulates the treatment sequence response probabilities for each treatment sequence as well as for the end of stage-1

library(LaplacesDemon)
PosteriorTrtSeqProbDesign1 <- function(niter, dat) {
  
  #Arguments
  #niter: number of draws from the posteriors of response probabilities
  #dat: dataset
  
  #End of study binary response indicator
  y <- dat$y
  
  #Stage-1 treatment assignment indicator
  a1 <- dat$a1
  
  #Stage-2 treatment assignment indicator
  a2 <- dat$a2
  
  #End of stage-1 binary response indicator
  s <- dat$s
  
  results <- matrix(NA, nrow = niter, ncol = 8)
  
  colnames(results) <- c("p_1","p_2","p_3","p_4","p_5","p_6","s1","s2")
  
  p_1_results <- p_2_results <- p_3_results <- p_4_results <- p_5_results <- p_6_results <- rep(0.5,niter)
  
  
  s1_results <- s2_results <- rep(0.5,niter)
  
  
  #Simulate from each of the six treatment sequences, the probability of response at the end of the trial from the posterior
  p_1_results <- rbeta(niter, shape1 = sum(y[a1==1&s==1])+1, sum(a1==1&s==1)-sum(y[a1==1&s==1])+1)
  p_2_results <- rbeta(niter, shape1 = sum(y[a1==1&s==0&a2==1])+1, sum(a1==1&s==0&a2==1)-sum(y[a1==1&s==0&a2==1])+1)
  p_3_results <- rbeta(niter, shape1 = sum(y[a1==1&s==0&a2==-1])+1, sum(a1==1&s==0&a2==-1)-sum(y[a1==1&s==0&a2==-1])+1)
  p_4_results <- rbeta(niter, shape1 = sum(y[a1==-1&s==1])+1, sum(a1==-1&s==1)-sum(y[a1==-1&s==1])+1)
  p_5_results <- rbeta(niter, shape1 = sum(y[a1==-1&s==0&a2==1])+1, sum(a1==-1&s==0&a2==1)-sum(y[a1==-1&s==0&a2==1])+1)
  p_6_results <- rbeta(niter, shape1 = sum(y[a1==-1&s==0&a2==-1])+1, sum(a1==-1&s==0&a2==-1)-sum(y[a1==-1&s==0&a2==-1])+1)
  
  #Simulate from each of the two first stage treatments the probability of response at the end of stage-1 from the posterior
  s1_results <- rbeta(niter, sum(s[a1==1])+1, sum(a1==1)-sum(s[a1==1])+1)
  s2_results <- rbeta(niter, sum(s[a1==-1])+1, sum(a1==-1)-sum(s[a1==-1])+1)
  
  
  results[,"p_1"] <- p_1_results
  results[,"p_2"]<-p_2_results
  results[,"p_3"] <- p_3_results
  results[,"p_4"]<-p_4_results
  results[,"p_5"]<-p_5_results
  results[,"p_6"]<-p_6_results
  results[,"s1"] <- s1_results
  results[,"s2"] <- s2_results
  
  return(results)
}
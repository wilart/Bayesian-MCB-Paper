#Compute the response probabilities for each treatment sequence for a given dataset

PosteriorTrtSeqProbGeneral <- function(niter, dat) {
  
  #Arguments
  #niter: number of draws from each posterior probability of response
  #dat: dataset
  
  
  results <- matrix(NA, nrow = niter, ncol = 10)
  
  colnames(results) <- c("p_1","p_2","p_3","p_4","p_5","p_6","p_7","p_8","s1","s2")
  
  #Final binary outcome
  y <- dat$y
  
  #Stage-1 treatment assignment
  a1<-dat$a1
  
  #Stage-2 treatment assignment for responders to stage-1 treatment
  a2r<-dat$a2r
  
  #Stage-2 treatment assignment for non-responders to stage-1 treatment
  a2nr<-dat$a2nr
  
  #Stage-1 response indicator (1 if responder at end of stage-1, 0 otherwise)
  s<-dat$s
  
  
  #Compute posterior probabilities of response at the end of the study for each of the eight embedded treatment sequences
  p_1_results <- rbeta(niter, shape1 = sum(y[a1==1&s==1&a2r==1]) + 1, shape2=sum(a1==1&s==1&a2r==1)-sum(y[a1==1&s==1&a2r==1])+1)
  p_2_results <- rbeta(niter, shape1 = sum(y[a1==1&s==1&a2r==-1]) + 1, shape2 = sum(a1==1&s==1&a2r==-1)-sum(y[a1==1&s==1&a2r==-1])+1)
  p_3_results <- rbeta(niter, shape1 = sum(y[a1==1&s==0&a2nr==1]) + 1, shape2 = sum(a1==1&s==0&a2nr==1)-sum(y[a1==1&s==0&a2nr==1])+1)
  p_4_results <- rbeta(niter, shape1 = sum(y[a1==1&s==0&a2nr==-1]) + 1, shape2 = sum(a1==1&s==0&a2nr==-1)-sum(y[a1==1&s==0&a2nr==-1])+1)
  p_5_results <- rbeta(niter, shape1 = sum(y[a1==-1&s==1&a2r==1]) + 1, shape2 = sum(a1==-1&s==1&a2r==1)-sum(y[a1==-1&s==1&a2r==1])+1)
  p_6_results <- rbeta(niter, shape1 = sum(y[a1==-1&s==1&a2r==-1]) + 1, shape2 = sum(a1==-1&s==1&a2r==-1)-sum(y[a1==-1&s==1&a2r==-1])+1)
  p_7_results <- rbeta(niter, shape1 = sum(y[a1==-1&s==0&a2nr==1]) + 1, shape2 = sum(a1==-1&s==0&a2nr==1)-sum(y[a1==-1&s==0&a2nr==1])+1)
  p_8_results <- rbeta(niter, shape1 = sum(y[a1==-1&s==0&a2nr==-1]) + 1, shape2 = sum(a1==-1&s==0&a2nr==-1)-sum(y[a1==-1&s==0&a2nr==-1])+1)
  
  #Compute posterior probability of response at the end of stage-1 for each of the two stage-1 treatments
  s1_results <- rbeta(niter, sum(s[a1==1]==1)+1, length(which(a1==1))-sum(s[a1==1])+1)
  s2_results <- rbeta(niter, sum(s[a1==-1]==1)+1, length(which(a1==-1))-sum(s[a1==-1])+1)
  
  
  results[,"p_1"] <- p_1_results
  results[,"p_2"]<-p_2_results
  results[,"p_3"] <- p_3_results
  results[,"p_4"]<-p_4_results
  results[,"p_5"]<-p_5_results
  results[,"p_6"]<-p_6_results
  results[,"p_7"] <- p_7_results
  results[,"p_8"] <- p_8_results
  
  results[,"s1"] <- s1_results
  results[,"s2"] <- s2_results
  return(results)
}


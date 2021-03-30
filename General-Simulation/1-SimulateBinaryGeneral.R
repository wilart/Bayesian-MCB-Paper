#Simulate n_sim data sets each with sample size given by "sample_size" for the general SMART study design

GeneralSimulateBinary <- function(sample_size=1000,
                                  n_sim,
                                  response_prob = c(0.5,0.9,0.7,0.6,0.8,0.8,0.4,0.5),
                                  stage_one_trt_one_response_prob = 0.7,
                                  stage_one_trt_two_response_prob = 0.4) {
  
  # Arguments
  # sample_size: sample size
  # n_sim: number of data set replicates
  # response_prob: vector of probabilities of response at the end of the study to each of the eight embedded treatment sequences
  # stage_one_trt_one_response_prob: probability of response to stage-1 treatment one
  # stage_one_trt_two_response_prob: probability of response to stage-1 treatment two
  
  sim_list <- vector("list", length = n_sim)
  
  for (i in 1:n_sim){

    #Stage-1 treatment assignment
    a1<-2*rbinom(sample_size,1,.5)-1 
    
    
    #Stage-1 responder and non-responder stage 2 treatment assignment
    a2r<-2*rbinom(sample_size,1,0.5)-1
    a2nr <- 2*rbinom(sample_size,1,0.5)-1
    
    #Response probabilities for stage-1
    s<-rep(NA,sample_size)
    s[a1==-1] <- rbinom(length(which(a1==-1)),size=1,stage_one_trt_one_response_prob)
    s[a1==1] <- rbinom(length(which(a1==1)),size=1,stage_one_trt_two_response_prob)
    #####################

    #Simulate binary outcomes
    y <- rep(NA,sample_size)
    
    y[a1==1&s==1&a2r==1] <- rbinom(length(which(a1==1&s==1&a2r==1)),1,prob = response_prob[1])
    y[a1==1&s==1&a2r==-1] <- rbinom(length(which(a1==1&s==1&a2r==-1)),1,prob = response_prob[2])
    y[a1==1&s==0&a2nr==1] <- rbinom(length(which(a1==1&s==0&a2nr==1)),1,prob = response_prob[3])
    y[a1==1&s==0&a2nr==-1] <- rbinom(length(which(a1==1&s==0&a2nr==-1)),1,prob = response_prob[4])
    y[a1==-1&s==1&a2r==1] <- rbinom(length(which(a1==-1&s==1&a2r==1)),1,prob = response_prob[5])
    y[a1==-1&s==1&a2r==-1] <- rbinom(length(which(a1==-1&s==1&a2r==-1)),1,prob = response_prob[6])
    y[a1==-1&s==0&a2nr==1] <- rbinom(length(which(a1==-1&s==0&a2nr==1)),1,prob = response_prob[7])
    y[a1==-1&s==0&a2nr==-1] <- rbinom(length(which(a1==-1&s==0&a2nr==-1)),1,prob = response_prob[8])
    
    
    sim_list[[i]]<-data.frame(a1,s,a2r,a2nr,y)
  }
  sim_list
}

set.seed(12374)

sim_binary_grid <- lapply(seq(150,500,50),function(x) GeneralSimulateBinary(sample_size=x, 
                                                                            n_sim=1000,
                                                                            response_prob = c(0.5,0.9,0.7,0.2,0.2,0.8,0.2,0.7),
                                                                            0.7,
                                                                            0.5))

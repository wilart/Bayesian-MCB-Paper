#Simulates n_sim datasets of size "sample_size" for a SMART design similar to an ENGAGE-type study. In the code, we call this 'Design-1'.

SimulateBinaryDesign1 <- 
  function(sample_size=1000,
           n_sim,
           response_prob=c(0.5,0.5,0.5,0.8,0.7,0.5),
           stage_one_trt_one_response_prob = 0.7,
           stage_one_trt_two_response_prob = 0.4){
  
    # Arguments:
    # sample_size: sample size
    # n_sim: number of replicate data sets
    # response_prob: vector of probabilities of response for each of the 6 embedded treatment sequences
    # stage_one_trt_one_response_prob: probability of response to first stage-1 treatment
    # stage_one_trt_two_response_prob: probability of response to second stage-1 treatment
    
    
  sim_list <- vector("list", length = n_sim)
  
  for (i in 1:n_sim){
    
    #First stage treatment indicator: coded as -1 and +1
    a1 <- 2*rbinom(sample_size,1,.5)-1 
    
    #Second stage treatment indicator: coded as -1 and +1
    a2 <- 2*rbinom(sample_size,1,0.5)-1
    
    #Stage-1 response probabilities
    s<-rep(NA,sample_size)
    s[a1==1] <- rbinom(length(which(a1==1)),size=1,stage_one_trt_one_response_prob)
    s[a1==-1] <- rbinom(length(which(a1==-1)),size=1,stage_one_trt_two_response_prob)

    #End-of-study outcomes
    y <- rep(NA,sample_size)
    y[a1==1&s==1] <- rbinom(length(which(a1==1&s==1)), size = 1, prob = response_prob[1])
    y[a1==1&s==0&a2==1] <- rbinom(length(which(a1==1&s==0&a2==1)), size = 1, prob = response_prob[2])
    y[a1==1&s==0&a2==-1] <- rbinom(length(which(a1==1&s==0&a2==-1)), size = 1, prob = response_prob[3])
    y[a1==-1&s==1] <- rbinom(length(which(a1==-1&s==1)), size = 1, prob = response_prob[4])
    y[a1==-1&s==0&a2==1] <- rbinom(length(which(a1==-1&s==0&a2==1)), size = 1, prob = response_prob[5])
    y[a1==-1&s==0&a2==-1] <- rbinom(length(which(a1==-1&s==0&a2==-1)), size = 1, prob = response_prob[6])
    
    
    sim_list[[i]]<-data.frame(a1,s,a2,y)
  }
  sim_list
  }

set.seed(3643)
sim_binary_grid <- lapply(seq(150,500,50),function(x) SimulateBinaryDesign1(sample_size=x, 
                                                                            n_sim=2000,
                                                                            response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
                                                                            0.7,
                                                                            0.5))

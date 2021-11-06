lambda_plus <- 0.5
lambda_minus <- 0.7


set.seed(12374)

sim_binary_grid <- lapply(seq(150,500,50),function(x) GeneralSimulateBinary(sample_size=x, 
                                                                            n_sim=1000,
                                                                            response_prob = c(0.5,0.9,0.7,0.2,0.2,0.8,0.2,0.7),
                                                                            0.5,
                                                                            0.7))
LogOR(c(0.5,0.9,0.7,0.2,0.2,0.8,0.2,0.7),0.7,0.5)






response_prob <- data.frame(expand.grid(c(0.8,0.9,0.95),seq(0.6,0.8,0.1),seq(0.2,0.5,0.1),seq(100,500,100)))
colnames(response_prob) <- c("phi_2",
                             "phi_3",
                             "phi",
                             "sample_size")
library(reshape2)

response_prob_melt <- melt(response_prob,id.vars="sample_size")

#######


stage_one_trt_one_response_prob <- 0.5
stage_one_trt_two_response_prob <- 0.7


########

#Compute power given inputs
ComputePowerBayesianGeneral <- function(sample_size=500,
                                        response_prob = c(0.5,0.9,0.7,0.2,0.2,0.8,0.2,0.7),
                                        stage_one_trt_one_response_prob = 0.5,
                                        stage_one_trt_two_response_prob = 0.7,
                                        threshold,
                                        alpha = 0.05
){
  
  #Arguments
  #sample_size: total sample size in SMART study
  #response_prob: probability of response for each of eight embedded treatment sequences
  #stage_one_trt_one_response_prob: probability of response at the end of stage-1 to treatment a1=1
  #stage_one_trt_two_response_prob: probability of response at the end of stage-1 to treatment a1=0
  #threshold: threshold of log-OR for exclusion from the set of best
  #alpha: probability of excluding true best EDTR from set of best
  
  pb <- txtProgressBar(min=0,
                       max=1000,
                       initial = 0,
                       style=3)
  
  upper_limit <- array(NA, dim=c(1000,10,8))
  
  
  
  #Stage-1 response indicator (1 if responder, 0 otherwise)
  s <- rep(NA,sample_size)
  
  
  for (i in 1:1000) {
    
    #Stage-1 treatment indicator
    a1 <- rbinom(sample_size,1,0.5)
    s[a1==1] <- rbinom(length(which(a1==1)),1,stage_one_trt_one_response_prob)
    s[a1==0] <- rbinom(length(which(a1==0)),1,stage_one_trt_two_response_prob)
    
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
      #Draw 1000 draws from the posterior of the probability of response at the end of the study for each of the eight embedded treatment sequences.
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
        logORThreshold <- LogOR(response_prob,
                                stage_one_trt_one_response_prob,
                                stage_one_trt_two_response_prob)
        thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
        max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
        
        log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=1000,ncol=8)
        
        
        #rank_matrix <- apply(log_OR_matrix,2,rank,ties.method = 'random')
        rank_matrix <- apply(log_OR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
        
        rank_max <- apply(rank_matrix,1,max)
        
        
        new_dat <- apply(log_OR_matrix[,],2,sort)
        
        ranks_quantile <- floor(quantile(rank_max,1-alpha))
        
        upper_limit[i,j,] <-new_dat[ranks_quantile,]
      
      
      
      
    }
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  
  
    rejection_indices <- which(abs(logORThreshold)>threshold)

  if (length(rejection_indices)==1) {
    return(mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0))
    
  } else {
    return(mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod)))
  }
}


library(doRNG)
library(doParallel)

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()



tst_0.5_interaction_start <- Sys.time()

set.seed(1274)
power_sensitivity_General_parameter4 <-foreach(z=1:180) %dorng% ComputePowerBayesianGeneral(sample_size = response_prob[z,4],
                                                                                                                                                         as.vector(c(response_prob[z,3],response_prob[z,1],response_prob[z,2],as.vector(response_prob[z,3]),(rep(response_prob[z,3],4))),mode="numeric"),
                                                                                                                                                         #response_prob =response_prob[[2]][[2]],
                                                                                                                                                         stage_one_trt_one_response_prob = stage_one_trt_one_response_prob,
                                                                                                                                                         stage_one_trt_two_response_prob = stage_one_trt_two_response_prob,
                                                                                                                                                         threshold = 0.0001)

stopCluster(cl)
Sys.time()
tst_0.5_interaction_stop <- Sys.time()

LogOR(response_prob=c(mean(dat_for_analysis_4_outcome_list[[z]][which(dat_for_analysis_4$a1==1&dat_for_analysis_4$s==1)]),
                      mean(dat_for_analysis_4_outcome_list[[z]][which(dat_for_analysis_4$a1==1&dat_for_analysis_4$s==0&dat_for_analysis_4$a2==1)]),
                      mean(dat_for_analysis_4_outcome_list[[z]][which(dat_for_analysis_4$a1==1&dat_for_analysis_4$s==0&dat_for_analysis_4$a2==-1)]),
                      mean(dat_for_analysis_4_outcome_list[[z]][which(dat_for_analysis_4$a1==-1&dat_for_analysis_4$s==1)]),
                      mean(dat_for_analysis_4_outcome_list[[z]][which(dat_for_analysis_4$a1==-1&dat_for_analysis_4$s==0&dat_for_analysis_4$a2==1)]),
                      mean(dat_for_analysis_4_outcome_list[[z]][which(dat_for_analysis_4$a1==-1&dat_for_analysis_4$s==0&dat_for_analysis_4$a2==-1)])),stage_one_trt_one_response_prob=lambda_plus,
      stage_one_trt_two_response_prob=lambda_minus)


LogOR(
                           as.vector(c(as.vector(response_prob[z,1:2]),(rep(response_prob[z,3],4))),mode="numeric"),
                           #response_prob =response_prob[[2]][[2]],
                           stage_one_trt_one_response_prob = stage_one_trt_one_response_prob,
                           stage_one_trt_two_response_prob = stage_one_trt_two_response_prob)


power_sensitivity_ENGAGE_parameter3[[10]] <- foreach(sample_size = seq(200,500,50)) %dopar% computePowerBayesian(sample_size = sample_size,
                                                                  response_prob[[10]],
                                                                  stage_one_trt_one_response_prob = stage_one_trt_one_response_prob,
                                                                  stage_one_trt_two_response_prob = stage_one_trt_two_response_prob,
                                                                  threshold = 0.7,
                                                                  alpha = 0.05)

library(reshape2)

library(tidyverse)

response_prob_mat <- cbind(response_prob,unlist(power_sensitivity_General_parameter4))
colnames(response_prob_mat) <- c("phi_2","phi_3","phi","sample_size", "Power")
power_sensitivity_melt <- melt(response_prob_mat,id.var=c("sample_size","Power","phi_2","phi_3"))


library(ggplot2)
power_sensitivity_melt$phi_2 <- as.factor(power_sensitivity_melt$phi_2)
power_sensitivity_melt$phi_3 <- as.factor(power_sensitivity_melt$phi_3)
power_sensitivity_melt$value <- as.factor(power_sensitivity_melt$value)
library(latex2exp)

levels(power_sensitivity_melt$phi_2) <- c(TeX('$\\phi_2 = 0.8$'),
                                          TeX('$\\phi_2 = 0.9$'),
                                          TeX('$\\phi_2 = 0.95$'))
levels(power_sensitivity_melt$phi_3) = c(TeX('$\\phi_3 = 0.6$'),
                                         TeX('$\\phi_3 = 0.7$'),
                                         TeX('$\\phi_3 = 0.8$'))
levels(power_sensitivity_melt$value) <- c(TeX('$\\phi_{inf} = 0.2'),
                                          TeX('$\\phi_{inf} = 0.3'),
                                          TeX('$\\phi_{inf} = 0.4'),
                                          TeX('$\\phi_{inf} = 0.5'))
labels <- unname(TeX(c(('$\\phi_{inf} = 0.2'),
                       ('$\\phi_{inf} = 0.3'),
                       ('$\\phi_{inf} = 0.4'),
                       ('$\\phi_{inf} = 0.5'))))
#pdf('Sensitivity-analysis-binary-outcome-General-11-6-21.pdf',width = 17,height=12)

ggplot2::ggplot(power_sensitivity_melt ,aes(x=`sample_size`,y=Power,colour=as.factor(value),linetype=as.factor(value)))+scale_y_continuous(breaks=seq(0,1,0.1))+geom_point(size=2)+geom_line(size=1,) + facet_grid((phi_3)~(phi_2),labeller=label_parsed)+labs(color="Inferior Prob. of Response",linetype="Inferior Prob. of Response")+ ggtitle("Power vs. Sample Size Sensitivity Analysis Binary Outcome, General Simulation") + scale_x_continuous(breaks=seq(100,500,100)) + theme_bw()+ theme(plot.title = element_text(hjust = 0.5,face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                          axis.ticks.length = unit(0.5, "cm"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                          axis.text=element_text(size=16), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                          axis.title.y=element_text(size=16,face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                          axis.title.x=element_text(size=16,face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                          legend.key.width = unit(1.5,"cm"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                          strip.text = element_text(size=18),
                                                                                                                                                                                                                                                                                                                                                                                                                                                          legend.text = element_text(size=18))+xlab("n")+ylab("Power")+
  scale_linetype_manual(name="Inferior Prob. of Response",values=c("twodash","dashed","solid","dotted"),labels=labels)+scale_color_discrete(labels=labels)

#dev.off() 
#save(power_sensitivity_melt,file="power_sensitivity_melt.rda")
###########

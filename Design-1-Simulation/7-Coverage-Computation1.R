args(SimulateBinaryDesign1)
set.seed(3752)
dat <- SimulateBinaryDesign1(sample_size = 150,
                      n_sim = 5000)





library(tidyverse)

computeEDTRDrawsFromSim <- function(x,dat,niter=10000) {
#bayesian_thetas_100 <- lapply(sim_1,function(x) ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=1000, dat= x)))
#bayesian_outcomes_100 <- colMeans(do.call(rbind,lapply(bayesian_thetas_100, function(x) colMeans(x))))

#bayesian_SE_100 <- apply(do.call(rbind,lapply(bayesian_thetas_100, function(x) colMeans(x))),2,sd)


beta_probit_current <- MASS::mvrnorm(1,c(0,0,0),diag(1,length(c(0,0,0))))

s_continuous <- rnorm(nrow(dat[[x]]),0,1)

beta_probit_results <- matrix(NA, nrow = niter,ncol= 3)
sigma_results <- rep(NA,niter)
sigma_probit_current <- 1

dat[[x]] <- dat[[x]] %>% filter(s == 1) %>% rbind(dat[[x]] %>% filter(s == 1) %>% mutate(a2 =-a2), dat[[x]] %>% filter(s==0))


w <- (dat[[x]]$s==1)/(0.5) + 
  (dat[[x]]$s==0)/(0.25)
y <- dat[[x]][,4]
dat[[x]] <- as.matrix(dat[[x]])[,-c(4)]



for (j in 1:niter) {
  
  for (i in 1:length(y)) {
    if (y[i]==1) {
      s_continuous[i] <- truncnorm::rtruncnorm(1,0,Inf,dat[[x]][i,] %*% beta_probit_current,1)
    } else {
      s_continuous[i] <- truncnorm::rtruncnorm(1,-Inf,0,dat[[x]][i,] %*% beta_probit_current,1)
      
    }
  }
  
  beta_probit_current <- MASS::mvrnorm(1,solve(t(dat[[x]])%*%diag(w)%*%dat[[x]])%*%t(dat[[x]])%*%diag(w)%*%s_continuous, sigma_probit_current* solve(t(dat[[x]])%*%diag(w)%*%dat[[x]]))
  
  
  RSS <- t(s_continuous-dat[[x]] %*% beta_probit_current)%*%diag(w)%*%(s_continuous-dat[[x]] %*% beta_probit_current)
  
  sigma_probit_current <- LaplacesDemon::rinvchisq(1,length(y)-3,RSS/(length(y)-3))
  
  
  beta_probit_results[j,] <- t(beta_probit_current)
  
  sigma_results[j] <- sigma_probit_current
  print(j)

}
return(beta_probit_results)
}

results <- computeEDTRDrawsFromSim(1)
EDTR_mat <- matrix(NA,nrow=niter,ncol=4)


EDTR_mat[,1] <- pnorm(results[,1]+results[,2]+results[,3])
EDTR_mat[,2] <- pnorm(results[,1]+results[,2]-results[,3])
EDTR_mat[,3] <- pnorm(results[,1]-results[,2]+results[,3])
EDTR_mat[,4] <- pnorm(results[,1]-results[,2]-results[,3])

#Credible interval coverage
thetadraws_log_odds <- log(EDTR_mat/(1-EDTR_mat))

#Compute index of best EDTR
max_odds_ind <- which.max(colMeans(thetadraws_log_odds))

#Compute log-odds ratios between each EDTR and best
Log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=nrow(thetadraws),ncol=4)


mean(apply(Log_OR_matrix<=matrix(ComputeMCBUpperLimitsDesign1(EDTR_mat,0.05),nrow=niter,ncol=4,byrow = T),1,prod))

##########
library(doParallel)
expit <- function(x) exp(x)/(1+exp(x))

#Number of cores/threads
no_cores <- detectCores()-1

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()

set.seed(1274)
clusterSetRNGStream(cl, 123)

tst_0.5_interaction_start <- Sys.time()
#47 replicates, 10000 iterations, 5 mixture components
#Runs main effects or interaction depending on what computeBetas is. That is, whether you run 2-ENGAGE-Interaction-Fit-Full-Conditionals.R or
#2-ENGAGE-Main-Effects-Fit-Full-Conditionals.R before. The object title can be chosen appropriately
results_dat_sim_coverage_5000_rep_10000iter <- foreach(i =1:5000,.errorhandling = "pass",.packages=c("dplyr","tmvtnorm",
                                                                                             "truncnorm",
                                                                                             "LaplacesDemon",
                                                                                             "MASS",
                                                                                             "condMVNorm","geepack")) %dopar% computeEDTRDrawsFromSim(i,dat,niter=10000)

stopCluster(cl)
Sys.time()

tst_0.5_interaction_stop <- Sys.time()



#############

#load("results_dat_sim_coverage_5000_rep_10000iter.rda")


computeEDTRsFromMSM <- function(results) {
  EDTR_mat <- matrix(NA,nrow=nrow(results),ncol=4)
  
EDTR_mat[,1] <- pnorm(results[,1]+results[,2]+results[,3])
EDTR_mat[,2] <- pnorm(results[,1]+results[,2]-results[,3])
EDTR_mat[,3] <- pnorm(results[,1]-results[,2]+results[,3])
EDTR_mat[,4] <- pnorm(results[,1]-results[,2]-results[,3])
return(EDTR_mat)
}

computeDifference <- function(thetadraws) {
  upper_limit <- rep(NA,4)
  thetadraws <- thetadraws
  #Compute log-OR
  thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
  
  #Compute index of best EDTR
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
  
  #Compute log-odds ratios between each EDTR and best
  Log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=nrow(thetadraws),ncol=4)
  return(Log_OR_matrix)
}


########
#Frequentist
differences  <- apply(do.call(rbind,lapply(results_dat_sim_coverage_5000_rep_10000iter,function(z) apply((computeDifference(computeEDTRsFromMSM(z[1:10000,]))),2,mean))),2,mean)


#Frequentist
c("Frequentist",round(mean(unlist(lapply(1:5000,function(x) prod((ComputeMCBUpperLimitsDesign1(computeEDTRsFromMSM(results_dat_sim_coverage_5000_rep_10000iter[[x]][1:10000,]))>=differences))))),4))


#Bayesian
c("Bayesian", round(median(unlist(lapply(1:5000,function(x) mean(apply(matrix((ComputeMCBUpperLimitsDesign1(computeEDTRsFromMSM(results_dat_sim_coverage_5000_rep_10000iter[[x]][1:10000,]))),nrow=10000,ncol=4,byrow=T)>=computeDifference(computeEDTRsFromMSM(results_dat_sim_coverage_5000_rep_10000iter[[x]][1:10000,])),1,prod))))),4))

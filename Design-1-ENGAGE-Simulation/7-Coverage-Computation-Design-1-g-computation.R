

set.seed(3752)
dat_200 <- SimulateBinaryDesign1(sample_size = 200,
                             n_sim = 1000)









set.seed(4763)
results_dat_sim_coverage_200_1000_rep_5000iter2 <- lapply(dat_200,function(x)(ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=5000, dat= x))))


computeDifference <- function(thetadraws) {
  upper_limit <- rep(NA,4)
  thetadraws <- thetadraws
  #Compute log-OR
  thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
  
  #Compute index of best EDTR
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
  
  #Compute log-odds ratios between each EDTR and best
  Log_OR_matrix <- thetadraws_log_odds[,]-matrix(thetadraws_log_odds[,max_odds_ind],nrow=nrow(thetadraws_log_odds),ncol=4)
  return(Log_OR_matrix)
}


########
#Frequentist
differences  <- apply(do.call(rbind,lapply(results_dat_sim_coverage_200_1000_rep_5000iter2,function(z) apply((computeDifference((z[seq(1,5000,1),]))),2,mean))),2,mean)


#Frequentist
c("Frequentist",round(mean(unlist(lapply(1:1000,function(x) prod((ComputeMCBUpperLimitsDesign1((results_dat_sim_coverage_200_1000_rep_5000iter2[[x]][seq(1,5000,1),]))>=differences))))),4))


#Bayesian
c("Bayesian", round(median(unlist(lapply(1:1000,function(x) mean(apply(matrix((ComputeMCBUpperLimitsDesign1((results_dat_sim_coverage_200_1000_rep_5000iter2[[x]][seq(1,5000,1),]))),nrow=5000,ncol=4,byrow=T)>=computeDifference((results_dat_sim_coverage_200_1000_rep_5000iter2[[x]][seq(1,5000,1),])),1,prod))))),4))




set.seed(37352)
dat_400 <- SimulateBinaryDesign1(sample_size = 400,
                             n_sim = 1000)



set.seed(4763)
results_dat_sim_coverage_400_1000_rep_5000iter2 <- lapply(dat_400,function(x)(ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=5000, dat= x))))



#Frequentist
differences  <- apply(do.call(rbind,lapply(results_dat_sim_coverage_400_1000_rep_5000iter2,function(z) apply((computeDifference((z[seq(1,5000,1),]))),2,mean))),2,mean)


#Frequentist
c("Frequentist",round(mean(unlist(lapply(1:1000,function(x) prod((ComputeMCBUpperLimitsDesign1((results_dat_sim_coverage_400_1000_rep_5000iter2[[x]][seq(1,5000,1),]))>=differences))))),4))


#Bayesian
c("Bayesian", round(median(unlist(lapply(1:1000,function(x) mean(apply(matrix((ComputeMCBUpperLimitsDesign1((results_dat_sim_coverage_400_1000_rep_5000iter2[[x]][seq(1,5000,1),]))),nrow=5000,ncol=4,byrow=T)>=computeDifference((results_dat_sim_coverage_400_1000_rep_5000iter2[[x]][seq(1,5000,1),])),1,prod))))),4))

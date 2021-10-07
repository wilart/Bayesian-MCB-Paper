
set.seed(37521)
dat_200 <- GeneralSimulateBinary(sample_size = 200,
                             n_sim = 2500)



results_dat_sim_general_coverage_200_2500_rep_5000iter <- lapply(dat_200,function(x)(ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=5000, dat= x))))

computeDifference <- function(thetadraws) {
  upper_limit <- rep(NA,8)
  thetadraws <- thetadraws
  #Compute log-OR
  thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
  
  #Compute index of best EDTR
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
  
  #Compute log-odds ratios between each EDTR and best
  Log_OR_matrix <- thetadraws_log_odds[,]-matrix(thetadraws_log_odds[,max_odds_ind],nrow=nrow(thetadraws_log_odds),ncol=8)
  
  return(Log_OR_matrix)
}


########
#Frequentist
differences  <- apply(do.call(rbind,lapply(results_dat_sim_general_coverage_200_2500_rep_5000iter,function(z) apply((computeDifference((z[seq(1,5000,1),]))),2,median))),2,median)


#Frequentist
c("Frequentist",round(mean(unlist(lapply(1:2500,function(x) prod((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_coverage_200_2500_rep_5000iter[[x]][seq(1,5000,1),]))>=differences))))),4))


#Bayesian
c("Bayesian", round(mean(unlist(lapply(1:2500,function(x) mean(apply(matrix((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_coverage_200_2500_rep_5000iter[[x]][seq(1,5000,1),]))),nrow=5000,ncol=8,byrow=T)>=computeDifference((results_dat_sim_general_coverage_200_2500_rep_5000iter[[x]][seq(1,5000,1),])),1,prod))))),4))




set.seed(373521)
dat_400 <- GeneralSimulateBinary(sample_size = 400,
                             n_sim = 2500)


results_dat_sim_general_coverage_400_2500_rep_5000iter <- lapply(dat_400,function(x)(ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=5000, dat= x))))

########
#Frequentist
differences  <- apply(do.call(rbind,lapply(results_dat_sim_general_coverage_400_2500_rep_5000iter,function(z) apply((computeDifference((z[seq(1,5000,1),]))),2,median))),2,median)


#Frequentist
c("Frequentist",round(mean(unlist(lapply(1:2500,function(x) prod((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_coverage_400_2500_rep_5000iter[[x]][seq(1,5000,1),]))>=differences))))),4))


#Bayesian
c("Bayesian", round(mean(unlist(lapply(1:2500,function(x) mean(apply(matrix((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_coverage_400_2500_rep_5000iter[[x]][seq(1,5000,1),]))),nrow=5000,ncol=8,byrow=T)>=computeDifference((results_dat_sim_general_coverage_400_2500_rep_5000iter[[x]][seq(1,5000,1),])),1,prod))))),4))

library(dplyr)

true_outcomes <- c(0.5*0.5+0.7*0.5,
                   0.5*0.5+0.2*0.5,
                   0.9*0.5+0.7*0.5,
                   0.9*0.5+0.2*0.5,
                   0.2*0.7+0.2*0.3,
                   0.2*0.7+0.7*0.3,
                   0.8*0.7+0.2*0.3,
                   0.8*0.7+0.7*0.3)



set.seed(3751)
dat_200 <- GeneralSimulateBinary(sample_size = 200,
                             n_sim = 400)

results_dat_sim_estimation_200_400_rep_5000iter_general <- lapply(dat_200,function(x)(ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=5000, dat= x))))

set.seed(37511)
dat_400 <- GeneralSimulateBinary(sample_size = 400,
                             n_sim = 400)



results_dat_sim_general_coverage_400_1000_rep_5000iter <- lapply(dat_400,function(x)(ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=5000, dat= x))))





MSM_outcomes_400 <- colMeans(do.call(rbind,lapply(results_dat_sim_general_coverage_400_1000_rep_5000iter,function(x)colMeans((x)))))
MSM_SE_400 <- apply(do.call(rbind,lapply(results_dat_sim_general_coverage_400_1000_rep_5000iter,function(x)colMeans((x)))),2,sd)
MSM_SD_400 <- apply(do.call(rbind,lapply(results_dat_sim_general_coverage_400_1000_rep_5000iter,function(x)apply((x),2,sd))),2,mean)


#######


MSM_outcomes_200 <- colMeans(do.call(rbind,lapply(results_dat_sim_estimation_200_400_rep_5000iter_general,function(x)colMeans((x)))))
MSM_SE_200 <- apply(do.call(rbind,lapply(results_dat_sim_estimation_200_400_rep_5000iter_general,function(x)colMeans((x)))),2,sd)
MSM_SD_200 <- apply(do.call(rbind,lapply(results_dat_sim_estimation_200_400_rep_5000iter_general,function(x)apply((x),2,sd))),2,mean)


knitr::kable(round(data.frame("Bias 200" = abs(MSM_outcomes_200-true_outcomes),
                              "MCSE 200" = MSM_SE_200,
                              "SD 200" = MSM_SD_200,
                              "Bias 400" = abs(MSM_outcomes_400-true_outcomes),
                              "MCSE 400" = MSM_SE_400,
                              "SD 400" =MSM_SD_400),4),format="latex",caption = "Simulation study")

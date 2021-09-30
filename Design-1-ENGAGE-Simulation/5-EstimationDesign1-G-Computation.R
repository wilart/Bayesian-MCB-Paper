

true_outcomes <- c(0.5*0.7+0.5*0.3,
                   0.5*0.7+0.5*0.3,
                   0.8*0.4+0.7*0.6,
                   0.8*0.4+0.5*0.6)
############
############

set.seed(37511)
dat_200 <- SimulateBinaryDesign1(sample_size = 200,
                             n_sim = 400)

set.seed(4763)
results_dat_sim_estimation_200_400_rep_5000iter_design_1 <- lapply(dat_200,function(x)(ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=5000, dat= x))))

set.seed(375112)
dat_400 <- SimulateBinaryDesign1(sample_size = 400,
                             n_sim = 400)

results_dat_sim_estimation_400_400_rep_5000iter_design_1 <- lapply(dat_400,function(x)(ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=5000, dat= x))))

set.seed(8735)
sim_true <- SimulateBinaryDesign1(sample_size=10000,n_sim = 10)




bayesian_thetas_400 <- colMeans(do.call(rbind,lapply(sim_true,function(x)colMeans(ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=5000, dat= x))))))


MSM_outcomes_400 <- colMeans(do.call(rbind,lapply(results_dat_sim_estimation_400_400_rep_5000iter_design_1,function(x)colMeans((x)))))
MSM_SE_400 <- apply(do.call(rbind,lapply(results_dat_sim_estimation_400_400_rep_5000iter_design_1,function(x)colMeans((x)))),2,sd)
MSM_SD_400 <- apply(do.call(rbind,lapply(results_dat_sim_estimation_400_400_rep_5000iter_design_1,function(x)apply((x),2,sd))),2,mean)


#######


MSM_outcomes_200 <- colMeans(do.call(rbind,lapply(results_dat_sim_estimation_200_400_rep_5000iter_design_1,function(x)colMeans((x)))))
MSM_SE_200 <- apply(do.call(rbind,lapply(results_dat_sim_estimation_200_400_rep_5000iter_design_1,function(x)colMeans((x)))),2,sd)
MSM_SD_200 <- apply(do.call(rbind,lapply(results_dat_sim_estimation_200_400_rep_5000iter_design_1,function(x)apply((x),2,sd))),2,mean)


knitr::kable(round(data.frame("Bias 200" = abs(MSM_outcomes_200-true_outcomes),
                              "MCSE 200" = MSM_SE_200,
                              "SD 200" = MSM_SD_200,
                              "Bias 400" = abs(MSM_outcomes_400-true_outcomes),
                              "MCSE 400" = MSM_SE_400,
                              "SD 400" =MSM_SD_400),4),format="latex",caption = "Simulation study")

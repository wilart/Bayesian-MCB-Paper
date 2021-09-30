





set.seed(37521)
dat_200 <- GeneralSimulateBinary(sample_size = 200,response_prob = rep(0.5,8), stage_one_trt_one_response_prob =0.5,stage_one_trt_two_response_prob = 0.5,
                             n_sim = 1000)


results_dat_sim_general_inclusion_200_5000_rep_5000iter <- lapply(dat_200,function(x)(ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=5000, dat= x))))


set.seed(375234)
dat_400 <- GeneralSimulateBinary(sample_size = 400,response_prob = rep(0.5,8), stage_one_trt_one_response_prob =0.5,stage_one_trt_two_response_prob = 0.5,
                             n_sim =1000)
##########
results_dat_sim_general_inclusion_400_1000_rep_5000iter <- lapply(dat_400,function(x)(ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=5000, dat= x))))



#Frequentist# Each EDTR is included in the set of best with probability 95% under null hypothesis
c("EDTR 1",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[1]>=0))),2))
c("EDTR 2",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[2]>=0))),2))
c("EDTR 3",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[3]>=0))),2))
c("EDTR 4",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[4]>=0))),2))
c("EDTR 5",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[5]>=0))),2))
c("EDTR 6",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[6]>=0))),2))
c("EDTR 7",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[7]>=0))),2))
c("EDTR 8",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_200_5000_rep_5000iter[[x]][seq(1,5000,1),]))))[8]>=0))),2))



c("EDTR 1",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[1]>=0))),2))
c("EDTR 2",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[2]>=0))),2))
c("EDTR 3",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[3]>=0))),2))
c("EDTR 4",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[4]>=0))),2))
c("EDTR 5",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[5]>=0))),2))
c("EDTR 6",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[6]>=0))),2))
c("EDTR 7",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[7]>=0))),2))
c("EDTR 8",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsGeneral((results_dat_sim_general_inclusion_400_1000_rep_5000iter[[x]][seq(1,5000,1),]))))[8]>=0))),2))


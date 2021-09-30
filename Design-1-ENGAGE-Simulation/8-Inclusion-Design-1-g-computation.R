




set.seed(3752)
dat_200 <- SimulateBinaryDesign1(sample_size = 200,response_prob = rep(0.5,6), stage_one_trt_one_response_prob =0.5,stage_one_trt_two_response_prob = 0.5,
                             n_sim = 1000)


set.seed(47631)
results_dat_sim_inclusion_200_1000_rep_5000iter2 <- lapply(dat_200,function(x)(ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=5000, dat= x))))



#######
set.seed(37521)
dat_400 <- SimulateBinaryDesign1(sample_size = 400,response_prob = rep(0.5,6), stage_one_trt_one_response_prob =0.5,stage_one_trt_two_response_prob = 0.5,
                                 n_sim = 1000)




set.seed(47631)
results_dat_sim_inclusion_400_1000_rep_5000iter2 <- lapply(dat_400,function(x)(ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=5000, dat= x))))


#############

##############


#Frequentist# Each EDTR is included in the set of best with probability approximately 95%
c("EDTR 1",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_200_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[1]>=0))),2))
c("EDTR 2",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_200_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[2]>=0))),2))
c("EDTR 3",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_200_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[3]>=0))),2))
c("EDTR 4",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_200_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[4]>=0))),2))



c("EDTR 1",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_400_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[1]>=0))),2))
c("EDTR 2",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_400_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[2]>=0))),2))
c("EDTR 3",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_400_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[3]>=0))),2))
c("EDTR 4",round(mean(unlist(lapply(1:1000,function(x) ((ComputeMCBUpperLimitsDesign1((results_dat_sim_inclusion_400_1000_rep_5000iter2[[x]][seq(1,5000,1),]))))[4]>=0))),2))


#Evaluate performance of Bayesian approach vs. marginal structural model.
#Produces the content of Table 3.
library(dplyr)

true_outcomes <- c(0.5*0.4+0.7*0.6,
                   0.5*0.4+0.6*0.6,
                   0.9*0.4+0.7*0.6,
                   0.9*0.4+0.6*0.6,
                   0.8*0.7+0.4*0.3,
                   0.8*0.7+0.5*0.3,
                   0.8*0.7+0.4*0.3,
                   0.8*0.7+0.5*0.3)

############
############

set.seed(2374)
sim_general <- GeneralSimulateBinary(sample_size=100,n_sim = 1000)

bayesian_thetas_100 <- lapply(sim_general,function(x) ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=10000, dat= x)))
bayesian_outcomes_100 <- colMeans(do.call(rbind,lapply(bayesian_thetas_100, function(x) colMeans(x))))

bayesian_SE_100 <- apply(do.call(rbind,lapply(bayesian_thetas_100, function(x) colMeans(x))),2,sd)

############################Misspecified MSM

cont <- rbind(c(1, 1, 1, 1),
              c(1, 1, 1,-1),
              c(1, 1,-1, 1),
              c(1, 1,-1,-1),
              c(1,-1, 1, 1),
              c(1,-1, 1,-1),
              c(1,-1,-1, 1),
              c(1,-1,-1,-1))
computeGLM <- function(dat) {
  dat_rep <- dat %>% filter(s == 1) %>% rbind(dat %>% filter(s == 1) %>% mutate(a2nr  = -a2nr),dat %>% filter(s==0) %>% rbind(dat%>%filter(s==0)%>%mutate(a2r=-a2r)))
 
  #w <- (dat_rep$s==1)/(1) +  
  #  (dat_rep$s==0&dat_rep$a1==-1)/(0.25)+(dat_rep$s==0&dat_rep$a1==1)
#  
  coef_fit <- coef(glm(y~a1+a2r+a2nr+a1*a2r+a1*a2nr,data=dat_rep,family=binomial))
  c(exp(coef_fit[1]+coef_fit[2]+coef_fit[3]+coef_fit[4]+coef_fit[5]+coef_fit[6])/(1+exp(coef_fit[1]+coef_fit[2]+coef_fit[3]+coef_fit[4]+coef_fit[5]+coef_fit[6])),
    exp(coef_fit[1]+coef_fit[2]+coef_fit[3]-coef_fit[4]+coef_fit[5]-coef_fit[6])/(1+exp(coef_fit[1]+coef_fit[2]+coef_fit[3]-coef_fit[4]+coef_fit[5]-coef_fit[6])),
    exp(coef_fit[1]+coef_fit[2]-coef_fit[3]+coef_fit[4]-coef_fit[5]+coef_fit[6])/(1+exp(coef_fit[1]+coef_fit[2]-coef_fit[3]+coef_fit[4]-coef_fit[5]+coef_fit[6])),
    exp(coef_fit[1]+coef_fit[2]-coef_fit[3]-coef_fit[4]-coef_fit[5]-coef_fit[6])/(1+exp(coef_fit[1]+coef_fit[2]-coef_fit[3]-coef_fit[4]-coef_fit[5]-coef_fit[6])),
    exp(coef_fit[1]-coef_fit[2]+coef_fit[3]+coef_fit[4]-coef_fit[5]-coef_fit[6])/(1+exp(coef_fit[1]-coef_fit[2]+coef_fit[3]+coef_fit[4]-coef_fit[5]-coef_fit[6])),
    exp(coef_fit[1]-coef_fit[2]+coef_fit[3]-coef_fit[4]-coef_fit[5]+coef_fit[6])/(1+exp(coef_fit[1]-coef_fit[2]+coef_fit[3]-coef_fit[4]-coef_fit[5]+coef_fit[6])),
    exp(coef_fit[1]-coef_fit[2]-coef_fit[3]+coef_fit[4]+coef_fit[5]-coef_fit[6])/(1+exp(coef_fit[1]-coef_fit[2]-coef_fit[3]+coef_fit[4]+coef_fit[5]-coef_fit[6])),
    exp(coef_fit[1]-coef_fit[2]-coef_fit[3]-coef_fit[4]+coef_fit[5]+coef_fit[6])/(1+exp(coef_fit[1]-coef_fit[2]-coef_fit[3]-coef_fit[4]+coef_fit[5]+coef_fit[6])))
}

MSM_outcomes_100 <- colMeans(do.call(rbind,lapply(sim_general, computeGLM)))
MSM_SE_100 <- apply(do.call(rbind,lapply(sim_general, computeGLM)),2,sd)
###############



set.seed(3957)
sim_general <- GeneralSimulateBinary(sample_size=400,n_sim = 1000)

bayesian_thetas_400 <- lapply(sim_general,function(x)ComputePosteriorEDTRProbsGeneral(PosteriorTrtSeqProbGeneral(niter=10000, dat= x)))


MSM_outcomes_400 <- colMeans(do.call(rbind,lapply(sim_general, computeGLM)))
MSM_SE_400 <- apply(do.call(rbind,lapply(sim_general, computeGLM)),2,sd)

bayesian_outcomes_400 <- colMeans(do.call(rbind,lapply(bayesian_thetas_400, function(x) colMeans(x))))

bayesian_SE_400 <- apply(do.call(rbind,lapply(bayesian_thetas_400, function(x) colMeans(x))),2,sd)


knitr::kable(round(data.frame("Bayesian Bias"=abs(bayesian_outcomes_100-true_outcomes),
                              "Bayesian SE"=bayesian_SE_100,
                              "MSM Bias"=abs(MSM_outcomes_100-true_outcomes),
                              MSM_SE_100,
                              "Bayesian Bias" = abs(bayesian_outcomes_400-true_outcomes),
                              bayesian_SE_400,
                              "MSM Bias" = abs(MSM_outcomes_400-true_outcomes),
                              MSM_SE_400),4),format="latex")

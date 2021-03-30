#Evaluate performance of the bayesian approach which uses Robin's G computation  procedure vs. a marginal structural model.
#Produces some of the content of Table 3

library(dplyr)

true_outcomes <- c(0.5*0.7+0.5*0.3,
                   0.5*0.7+0.5*0.3,
                   0.8*0.4+0.7*0.6,
                   0.8*0.4+0.5*0.6)
############
############

set.seed(2374)
sim_1 <- SimulateBinaryDesign1(sample_size=150,n_sim = 1000)


#load("results_MSM_niter_5000_1000_rep_400_n2.rda")
#############
set.seed(3826)
estimateMSMProbit <- function(dat,niter=5000) {

beta_probit_current <- MASS::mvrnorm(1,c(0,0,0,0),diag(1,length(c(0,0,0,0))))

s_continuous <- rnorm(nrow(dat),0,1)

beta_probit_results <- matrix(NA, nrow = niter,ncol= 4)
X_design_matrix <- data.frame(dat)
X_design_matrix <- X_design_matrix %>% filter(s == 1) %>% rbind(X_design_matrix %>% filter(s == 1) %>% mutate(a2 =-a2), X_design_matrix %>% filter(s==0))%>% mutate("a1_by_a2" = a1*a2)


w <- (X_design_matrix$s==1)/(0.5) + 
  (X_design_matrix$s==0)/(0.25)
y <- X_design_matrix[,4]
X_design_matrix <- as.matrix(X_design_matrix)[,-c(4)]
X_design_matrix <- as.matrix(cbind(1,X_design_matrix[,-2]))
for (j in 1:niter) {
  
  for (i in 1:length(y)) {
    if (y[i]==1) {
      s_continuous[i] <- truncnorm::rtruncnorm(1,0,Inf,X_design_matrix[i,] %*% beta_probit_current,1)
    } else {
      s_continuous[i] <- truncnorm::rtruncnorm(1,-Inf,0,X_design_matrix[i,] %*% beta_probit_current,1)
      
    }
  }
  
  beta_probit_current <- MASS::mvrnorm(1,solve(t(X_design_matrix)%*%diag(w)%*%X_design_matrix)%*%t(X_design_matrix)%*%diag(w)%*%s_continuous, solve(t(X_design_matrix)%*%diag(w)%*%X_design_matrix))
  
  


  
  beta_probit_results[j,] <- t(beta_probit_current)
  
  print(j)
}
colnames(beta_probit_results) <- c("intercept","a1","a2","a1_by_a2")

return(beta_probit_results)
}


colMeans(beta_probit_results[1:niter,])
apply(beta_probit_results,2,median)



library(doParallel)

#Number of cores/threads
no_cores <- detectCores()-1

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()

set.seed(1274)
clusterSetRNGStream(cl, 123)

tst_0.5_interaction_start <- Sys.time()
#save(results_MSM_niter_5000_1000_rep_400_n3,file="results_MSM_niter_5000_1000_rep_400_n3.rda")
#47 replicates, 10000 iterations, 5 mixture components
#Runs main effects or interaction depending on what computeBetas is. That is, whether you run 2-ENGAGE-Interaction-Fit-Full-Conditionals.R or
#2-ENGAGE-Main-Effects-Fit-Full-Conditionals.R before. The object title can be chosen appropriately
results_MSM_niter_5000_1000_rep_150_n3 <- foreach(i =1:1000,.errorhandling = "pass",.packages=c("dplyr","tmvtnorm",
                                                                                                     "truncnorm",
                                                                                                     "LaplacesDemon",
                                                                                                     "MASS")) %dopar% estimateMSMProbit(sim_1[[i]],niter=5000)

stopCluster(cl)
Sys.time()

tst_0.5_interaction_stop <- Sys.time()
#save(results_MSM_niter_5000_1000_rep_150_n3,file="results_MSM_niter_5000_1000_rep_150_n3.rda")
#save(results_MSM_niter_5000_1000_rep_400_n,file="results_MSM_niter_5000_1000_rep_400_n.rda")
#load("results_MSM_niter_5000_1000_rep_400_n.rda")
############
#Check
dat <- dat %>% filter(s == 1) %>% rbind(dat %>% filter(s == 1) %>% mutate(a2 =-a2), dat %>% filter(s==0))

w2 <- (dat$s==1)/(0.5) + 
  (dat$s==0)/(0.25)

coef_fit <- coef(glm(y~a1+a2+a1*a2,data=dat,weights=w2,family=binomial(link="probit")))
############################ MSM
#save(results_MSM_niter_5000_1000_rep_150_n,file="results_MSM_niter_5000_1000_rep_150_n.rda")
###############


#set.seed(3957)
#sim_1 <- SimulateBinaryDesign1(sample_size=400,n_sim = 1000)

#bayesian_thetas_400 <- lapply(sim_1,function(x)ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=1000, dat= x)))


#MSM_outcomes_400 <- colMeans(do.call(rbind,lapply(sim_1, computeGLM)))
#MSM_SE_400 <- apply(do.call(rbind,lapply(sim_1, computeGLM)),2,sd)

#bayesian_outcomes_400 <- colMeans(do.call(rbind,lapply(bayesian_thetas_400, function(x) colMeans(x))))

#bayesian_SE_400 <- apply(do.call(rbind,lapply(bayesian_thetas_400, function(x) colMeans(x))),2,sd)


#knitr::kable(round(data.frame("Bayesian Bias"=abs(bayesian_outcomes_100-true_outcomes),
#      "Bayesian SE"=bayesian_SE_100,
#      "MSM Bias"=abs(MSM_outcomes_100-true_outcomes),
#      MSM_SE_100,
#      "Bayesian Bias" = abs(bayesian_outcomes_400-true_outcomes),
#      bayesian_SE_400,
#      "MSM Bias" = abs(MSM_outcomes_400-true_outcomes),
#      MSM_SE_400),4),format="latex")
#####################

EDTR_outcomes <- colMeans(pnorm(cbind(results_MSM_niter_5000_1000_rep_400_n[[1]]
%*%c(1,1,1,1),
results_MSM_niter_5000_1000_rep_400_n[[1]]
%*%c(1,1,-1,-1),
results_MSM_niter_5000_1000_rep_400_n[[1]]
%*%c(1,-1,1,-1),
results_MSM_niter_5000_1000_rep_400_n[[1]]
%*%c(1,-1,-1,1))))


EDTR_outcomes <- colMeans(pnorm(cbind(beta_probit_results
                                      %*%c(1,1,1,1),
                                      beta_probit_results
                                      %*%c(1,1,-1,-1),
                                      beta_probit_results
                                      %*%c(1,-1,1,-1),
                                      beta_probit_results
                                      %*%c(1,-1,-1,1))))

colnames(EDTR_outcomes) <- c('EDTR 1',"EDTR 2","EDTR 3","EDTR 4")
############################ MSM
computeGLM <- function(dat) {
  dat <- dat %>% filter(s == 1) %>% rbind(dat %>% filter(s == 1) %>% mutate(a2 =-a2), dat %>% filter(s==0))
  
  w <- (dat$s==1)/(0.5) + 
    (dat$s==0)/(0.25)
  
  coef_fit <- coef(glm(y~a1+a2+a1*a2,data=dat,weights=w,family=binomial))
  c(exp(coef_fit[1]+coef_fit[2]+coef_fit[3]+coef_fit[4])/(1+exp(coef_fit[1]+coef_fit[2]+coef_fit[3]+coef_fit[4])),
    exp(coef_fit[1]+coef_fit[2]-coef_fit[3]-coef_fit[4])/(1+exp(coef_fit[1]+coef_fit[2]-coef_fit[3]-coef_fit[4])),
    exp(coef_fit[1]-coef_fit[2]+coef_fit[3]-coef_fit[4])/(1+exp(coef_fit[1]-coef_fit[2]+coef_fit[3]-coef_fit[4])),
    exp(coef_fit[1]-coef_fit[2]-coef_fit[3]+coef_fit[4])/(1+exp(coef_fit[1]-coef_fit[2]-coef_fit[3]+coef_fit[4])))
}
MSM_outcomes_100 <- colMeans(do.call(rbind,lapply(sim_1, computeGLM)))
MSM_SE_100 <- apply(do.call(rbind,lapply(sim_1, computeGLM)),2,sd)
###############



set.seed(3957)
sim_1 <- SimulateBinaryDesign1(sample_size=400,n_sim = 1000)

set.seed(8735)
sim_true <- SimulateBinaryDesign1(sample_size=10000,n_sim = 10)




bayesian_thetas_400 <- lapply(sim_1,function(x)ComputePosteriorEDTRProbsDesign1(PosteriorTrtSeqProbDesign1(niter=1000, dat= x)))


MSM_outcomes_400 <- colMeans(do.call(rbind,lapply(sim_1, computeGLM)))
MSM_SE_400 <- apply(do.call(rbind,lapply(sim_1, computeGLM)),2,sd)

bayesian_outcomes_400 <- colMeans(do.call(rbind,lapply(bayesian_thetas_400, function(x) colMeans(x))))

bayesian_SE_400 <- apply(do.call(rbind,lapply(bayesian_thetas_400, function(x) colMeans(x))),2,sd)


colMeans(do.call(rbind,lapply(1:1000,function(x) pnorm(t(cbind(c(1,1,1,1),
                                c(1,1,-1,-1),
                                c(1,-1,1,-1),
                                c(1,-1,-1,1))%*%colMeans(results_MSM_niter_5000_1000_rep_400_n[[x]]))))))
load("results_MSM_niter_5000_1000_rep_400_n3.rda")
load("results_MSM_niter_5000_1000_rep_150_n3.rda")
options(scipen=999)
knitr::kable(round(data.frame(
      "MSM Bias 150" = abs(colMeans(do.call(rbind,lapply(1:1000,function(x) colMeans(pnorm((results_MSM_niter_5000_1000_rep_150_n3[[x]])%*%cbind(c(1,1,1,1),
                                                                                                                                                 c(1,1,-1,-1),
                                                                                                                                                 c(1,-1,1,-1),
                                                                                                                                                 c(1,-1,-1,1)))))))-true_outcomes),
      "MSM_SE_150" = apply(do.call(rbind,lapply(1:1000,function(x) colMeans(pnorm((results_MSM_niter_5000_1000_rep_150_n3[[x]])%*%cbind(c(1,1,1,1),
                                                                                                                                           c(1,1,-1,-1),
                                                                                                                                           c(1,-1,1,-1),
                                                                                                                                           c(1,-1,-1,1)))))),2,sd),
    
      "MSM Bias 400" = abs(colMeans(do.call(rbind,lapply(1:1000,function(x) colMeans(pnorm((results_MSM_niter_5000_1000_rep_400_n3[[x]])%*%cbind(c(1,1,1,1),
                                                                                                                                           c(1,1,-1,-1),
                                                                                                                                           c(1,-1,1,-1),
                                                                                                                                           c(1,-1,-1,1)))))))-true_outcomes),
      apply(do.call(rbind,lapply(1:1000,function(x) colMeans(pnorm((results_MSM_niter_5000_1000_rep_400_n3[[x]])%*%cbind(c(1,1,1,1),
                                                                                                                                c(1,1,-1,-1),
                                                                                                                                c(1,-1,1,-1),
                                                                                                                                c(1,-1,-1,1)))))),2,sd)),4),format="latex")

      
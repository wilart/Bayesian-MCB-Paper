#Evaluate power prediction performance for the general SMART

n_grid <- seq(150,500,50)
mc_param_sim <- MCBUpperLimits <-vector("list",length = length(n_grid))


library(doParallel)
library(doRNG)

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()

n <- 2000; p <- length(n_grid)
rng <- RNGseq( n * p, 123)
tst_0.5_interaction_start <- Sys.time()
mc_param_sim <- foreach(j=1:2000)%:% foreach(i =1:length(n_grid),r=rng[(j-1)*p + 1:p],.errorhandling = "pass",.packages=c("dplyr",
                                                                                        "LaplacesDemon",
                                                                                        "MASS")) %dopar% {# set RNG seed
                                                                                          rngtools::setRNG(r)
                                                                                          PosteriorTrtSeqProbGeneral(niter=1000, sim_binary_grid[[i]][[j]])
}
stopCluster(cl)
Sys.time()

tst_0.5_interaction_stop <- Sys.time()
##########
library(doParallel)

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()

n <- 2000; p <- length(n_grid)
rng <- RNGseq( n * p, 123)
tst_0.5_interaction_start <- Sys.time()
MCBUpperLimits <- foreach(j=1:2000)%:% foreach(i =1:length(n_grid),r=rng[(j-1)*p + 1:p],.errorhandling = "pass",.packages=c("dplyr",
                                                                                                     "LaplacesDemon",
                                                                                                     "MASS")) %dopar% {
                                                                                                       # set RNG seed
                                                                                                       rngtools::setRNG(r)
                                                                                                       ComputeMCBUpperLimitsGeneral(ComputePosteriorEDTRProbsGeneral(mc_param_sim[[j]][[i]]))
}
stopCluster(cl)
Sys.time()

tst_0.5_interaction_stop <- Sys.time()



library(doParallel)

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()

rng <- RNGseq(length(n_grid),123)

tst_0.5_interaction_start <- Sys.time()
predicted_power_grid <- foreach(i =1:length(n_grid),r=rng[1:length(n_grid)],.errorhandling = "pass",.packages=c("dplyr",
                                                                                                    "LaplacesDemon",
                                                                                                    "MASS")) %dopar% { 
                                                                                                      rngtools::setRNG(r)
                                                                                                      
                                                                                                      ComputePowerBayesianGeneral(n_grid[i],
                                                                                                                                                 response_prob = c(0.5,0.9,0.7,0.2,0.2,0.8,0.2,0.7),
                                                                                                                                                 0.5,
                                                                                                                                                 0.7,
                                                                                                                                                 threshold = 0.8)
}
stopCluster(cl)
Sys.time()

tst_0.5_interaction_stop <- Sys.time()


summary_melt <- reshape2::melt(data.frame('n_grid'=n_grid,'Predicted Power'=unlist(predicted_power_grid),
                     'Empirical Power' = unlist(lapply(1:length(n_grid),function(w) mean(apply(do.call(rbind,lapply(MCBUpperLimits,function(x) x[[w]]))[,c(1,2,4,5,6,7)],1,function(z)prod(z<=0)))))),id.vars='n_grid')



library(ggplot2)
g_power_general <- ggplot(summary_melt,aes(x=n_grid,y=value,color=factor(variable),linetype=factor(variable)))+
  geom_line(size=0.8)+
  geom_point(size=4)+
  geom_hline(yintercept=0.8)+
  scale_y_continuous(breaks=seq(0,1,0.1),expand = c(0,0),limits = c(0,1))+
  scale_x_continuous(breaks=seq(150,500,50))+
  theme_bw()+
  xlab("n")+
  ylab("Power")+
  ggtitle("Power vs. Sample Size for General SMART with 8 EDTRs")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.ticks.length = unit(0.5, "cm"),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14,face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        legend.text=element_text(size=14),
        legend.key.width = unit(2,"cm"),
        legend.title=element_blank(),
        legend.position="none")

#pdf("General-power-plot-9-29-21.pdf")
g_power_general
#dev.off()


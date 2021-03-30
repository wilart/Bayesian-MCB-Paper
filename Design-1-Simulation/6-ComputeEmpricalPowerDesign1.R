#Evaluate power prediction performance.
#Produces the right-hand plot of Figure 3.
n_grid <- seq(150,500,50)

mc_param_sim <- tst <-vector("list",length = length(n_grid))

#Take simulated SMART dataset, compute the marginal response probabilities, transform to EDTR response probabilities, compute empirical power
set.seed(3743)
for (i in 1:length(seq(150,500,50))) {

  for (j in 1:1000) {

    
    mc_param_sim[[i]][[j]] <- PosteriorTrtSeqProbDesign1(niter=1000, sim_binary_grid[[i]][[j]])
    tst[[i]][[j]] <- ComputeMCBUpperLimitsDesign1(ComputePosteriorEDTRProbsDesign1(mc_param_sim[[i]][[j]]))
    
    print(j)
    
  }
  
  
}

#Print the log odds ratio
-round(LogOR(),2)



predicted_power_grid <- rep(NA,length(n_grid))
set.seed(3754)
for (i in 1:length(n_grid))  {
predicted_power_grid[i] <- ComputePowerBayesianDesign1(n_grid[i],
                     response_prob = c(0.5,0.9,0.3,0.7,0.5,0.8),
                     0.7,
                     0.5)
}


#Plot power as a function of sample size, both predicted and empirical
library(ggplot2)
g_power_design_1 <- ggplot(data.frame(n_grid, reshape2::melt(cbind('Predicted Power'=predicted_power_grid,'Empirical Power' = unlist(lapply(1:length(n_grid),function(x) mean(do.call(rbind,tst[[x]])[,2]<=0&do.call(rbind,tst[[x]])[,3]<=0)))))),aes(x=n_grid,y=value,color=Var2,linetype=Var2))+
  geom_line(size=0.8)+
  geom_point(size=4)+
  geom_hline(yintercept=0.8)+
  scale_y_continuous(breaks=seq(0,1,0.1),expand = c(0,0),limits = c(0,1))+
  scale_x_continuous(breaks=seq(150,500,50))+
  theme_bw()+
  xlab("n")+
  ylab("Power")+
  ggtitle("Power vs. Sample Size for SMART Design 1 with 4 EDTRs")+
  theme(plot.title = element_text(hjust = 0.5,face="bold"),
        axis.ticks.length = unit(0.5, "cm"),
        axis.text=element_text(size=14), 
        axis.title.y=element_text(size=14,face="bold"),
        axis.title.x=element_text(size=14,face="bold"),
        legend.text=element_text(size=14),
        legend.key.width = unit(2,"cm"),
        legend.title=element_blank(),
        legend.position="none") + 
  guides(linetype = guide_legend(override.aes = list(size = 1.5)))

pdf("design-1-power-plot-1-2-1-21.pdf")
g_power_design_1
dev.off()

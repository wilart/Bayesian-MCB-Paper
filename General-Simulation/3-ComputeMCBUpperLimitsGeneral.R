#General SMART: 8 treatment sequences and 8 embedded dynamic treatment regime.
#Compute simultaneous alpha level credible intervals from probabilities of response for each embedded dynamic treatment regime


ComputeMCBUpperLimitsGeneral <- function(thetadraws, alpha = 0.05) {
  
  #Arguments
  #thetadraws: draws of embedded dynamic treatment regime response probabilities
  #alpha: type I error rate (excluding the optimal embedded dynamic treatment regime)
  
  
  upper_limit <- rep(NA,8)

  #Compute log-OR
  thetadraws_log_odds <- log(thetadraws/(1-thetadraws))
  
  #Compute index of best EDTR
  max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
  
  #Compute log-odds ratios between each EDTR and best
  Log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=nrow(thetadraws),ncol=8)

  #Rank log-OR
  rank_matrix <- apply(Log_OR_matrix,2,rank,ties.method = 'random')

  #Find max rank
  rank_max <- apply(rank_matrix,1,max)

  #Create sorted log-OR
  new_dat <- apply(Log_OR_matrix[,],2,sort)

  #Compute 100(1-alpha)% upper quantile
  ranks_quantile <- ceiling(quantile(rank_max,1-alpha))

  #Compute upper limit of credible interval. One for each log-OR which determines the set of best.
  upper_limit <-new_dat[ranks_quantile,]


  return(upper_limit)
}

#Compute Log-Odds Ratios for each Embedded Dynamic Treatment Regime

LogOR <- function(response_prob = c(0.5,0.9,0.7,0.2,0.2,0.8,0.2,0.7),
                  stage_one_trt_one_response_prob = 0.5,
                  stage_one_trt_two_response_prob = 0.7) {
  
  # Arguments:
  # response_prob: probability of response for each of embedded treatment sequences. 
  # stage_one_trt_one_response_prob: probability of response to stage-1 treatment given initial treatment is 1
  # stage_one_trt_two_response_prob: probability of response to stage-1 treatment given initial treatment is 0


    EDTRs <- c(
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[3] * (1 - stage_one_trt_one_response_prob),
      response_prob[1] * stage_one_trt_one_response_prob + response_prob[4] * (1 - stage_one_trt_one_response_prob),
      response_prob[2] * stage_one_trt_one_response_prob + response_prob[3] * (1 - stage_one_trt_one_response_prob),
      response_prob[2] * stage_one_trt_one_response_prob + response_prob[4] * (1 - stage_one_trt_one_response_prob),
      response_prob[5] * stage_one_trt_two_response_prob + response_prob[7] * (1 - stage_one_trt_two_response_prob),
      response_prob[5] * stage_one_trt_two_response_prob + response_prob[8] * (1 - stage_one_trt_two_response_prob),
      response_prob[6] * stage_one_trt_two_response_prob + response_prob[7] * (1 - stage_one_trt_two_response_prob),
      response_prob[6] * stage_one_trt_two_response_prob + response_prob[8] * (1 - stage_one_trt_two_response_prob)
    )



  # Compute log-OR
  thetadraws_log_odds <- log(EDTRs / (1 - EDTRs))

  # Compute index of best EDTR
  max_odds_ind <- (which.max((thetadraws_log_odds)))

  # Compute log-odds ratios between each EDTR and best
  Log_OR_output <- matrix((thetadraws_log_odds - thetadraws_log_odds[max_odds_ind]), nrow = 1, ncol = length(thetadraws_log_odds))


    colnames(Log_OR_output) <- c(
      "EDTR 1",
      "EDTR 2",
      "EDTR 3",
      "EDTR 4",
      "EDTR 5",
      "EDTR 6",
      "EDTR 7",
      "EDTR 8"
    )


  return(Log_OR_output)
}

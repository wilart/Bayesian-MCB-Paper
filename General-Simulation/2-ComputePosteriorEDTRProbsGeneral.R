#Computes the embedded dynamic treatment regime (EDTR) response probability draws using Robins' G-computation formula

ComputePosteriorEDTRProbsGeneral <- function(x) {
  cbind(x[,"p_1"]*x[,"s1"]+x[,"p_3"]*(1-x[,"s1"]),
        x[,"p_1"]*x[,"s1"]+x[,"p_4"]*(1-x[,"s1"]),
        x[,"p_2"]*x[,"s1"]+x[,"p_3"]*(1-x[,"s1"]),
        x[,"p_2"]*x[,"s1"]+x[,"p_4"]*(1-x[,"s1"]),
        x[,"p_5"]*x[,"s2"]+x[,"p_7"]*(1-x[,"s2"]),
        x[,"p_5"]*x[,"s2"]+x[,"p_8"]*(1-x[,"s2"]),
        x[,"p_6"]*x[,"s2"]+x[,"p_7"]*(1-x[,"s2"]),
        x[,"p_6"]*x[,"s2"]+x[,"p_8"]*(1-x[,"s2"]))
}
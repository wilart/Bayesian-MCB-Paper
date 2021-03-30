#Apply Robin's G-computation formula to compute the embedded dynamic treatment regime (EDTR) draws
ComputePosteriorEDTRProbsDesign1 <- function(x) {
  cbind(x[,"p_1"]*(x[,"s1"])+x[,"p_2"]*(1-(x[,"s1"])),
        x[,"p_1"]*(x[,"s1"])+x[,"p_3"]*(1-(x[,"s1"])),
        x[,"p_4"]*(x[,"s2"])+x[,"p_5"]*(1-(x[,"s2"])),
        x[,"p_4"]*(x[,"s2"])+x[,"p_6"]*(1-(x[,"s2"])))
}
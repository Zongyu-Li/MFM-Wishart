# The function is adapted from the code for Yin et al. (2023)
# See MFM_MN.R in https://github.com/fyin-stats/MFM-MxN
# Fan Yin, Guanyu Hu & Weining Shen (2023) Analysis of Professional
# Basketball Field Goal Attempts via a Bayesian Matrix Clustering Approach, Journal of
# Computational and Graphical Statistics, 32:1, 49-60, DOI: 10.1080/10618600.2022.2085727

getDahl <- function(fit, burn){
  ##############################################################################
  ## Input: fit  = result from MFM_Wishart_collapsed or DPM_Wishart_collapsed ##
  ##        burn = the number of burn-in iterations                           ##
  ## Output:                                                                  ##
  ##        zout  = estimated clustering configuration, a n by 1 vector       ##
  ##        nuout = estimated shared degrees of freedom, a scalar             ##
  ##############################################################################
  iters <- fit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

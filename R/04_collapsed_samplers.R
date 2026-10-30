# The code adapted from the code for Yin et al. (2023)
# See MFM_MN.R in https://github.com/fyin-stats/MFM-MxN
# Fan Yin, Guanyu Hu & Weining Shen (2023) Analysis of Professional
# Basketball Field Goal Attempts via a Bayesian Matrix Clustering Approach, Journal of
# Computational and Graphical Statistics, 32:1, 49-60, DOI: 10.1080/10618600.2022.2085727

MFM_Wishart_collapsed <- function(
    data, niterations,
    GAMMA,            # Dirichlet gamma
    Psi0, kappa0,     # IW prior hyperparams
    nu_init, nu_L, nu_U,  # nu prior bounds + init
    initNClusters,    # Number of initial clusters
    proposal_sd = 1,  # RW-MH sd for each cluster-specific nu_c
    lambda = 1,       # K - 1 ~ Poisson(lambda)
    k_max = 500,      # truncation for VN
    m_aux = 1         # number of Neal Algorithm 8 auxiliary nu values
){
  
  ## data: p x p x n array, each slice SPD
  n <- dim(data)[3]
  p <- nrow(data[,,1])
  
  if(m_aux < 1 || m_aux != as.integer(m_aux)){
    stop("m_aux must be a positive integer.")
  }
  m_aux <- as.integer(m_aux)
  
  ## Precompute logdet(Wi)(for updating z and nu)
  logdetW <- rep(NA_real_, n)
  for(i in 1:n){
    logdetW[i] <- logdet_spd(data[,,i])
  }
  
  ## Precompute VN in the MFM prior (k truncated at a large integer)
  N <- n
  VN <- 0
  tmax <- min(n + 10, k_max)
  VN <- rep(NA_real_, tmax)
  
  for (t in 1:tmax){
    r <- log(0)
    for (k in t:k_max){
      b <- sum(log((k-t+1):k)) -
        sum(log((k*GAMMA):(k*GAMMA+N-1))) +
        dpois(k-1, lambda, log = TRUE)
      m <- max(b, r)
      r <- log(exp(r-m) + exp(b-m)) + m
    }
    VN[t] <- r
  }
  
  ## Initialize clusterAssign
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  ## for collapsed sampling
  ## For fast computation, we keep the sufficient stats in each cluster
  ##                S_c = sum W_i,  A_c = sum log|W_i|               ##
  clusterSizes <- table(as.factor(clusterAssign))
  nClusters <- length(clusterSizes)
  
  Sarr <- array(0, dim=c(p,p,nClusters))
  A    <- rep(0, nClusters)
  for(i in 1:n){
    c <- clusterAssign[i]
    Sarr[,,c] <- Sarr[,,c] + data[,,i]
    A[c] <- A[c] + logdetW[i]
  }
  
  ## A scalar nu_init is replicated; alternatively pass a vector of length initNClusters.
  if(length(nu_init) == 1){
    nu <- rep(nu_init, nClusters)
  } else {
    if(length(nu_init) != nClusters){
      stop("nu_init must be a scalar or a vector of length initNClusters.")
    }
    nu <- nu_init
  }
  
  History <- vector("list", niterations)
  
  ## start MCMC
  for(iter in 1:niterations){
    
    #if(iter %% 1000 == 0){
    #  print(iter)
    #}
    
    ## update z, the labels for each sample
    clusterSizes <- table(as.factor(clusterAssign))
    nClusters <- length(clusterSizes)
    
    for(i in 1:n){
      
      cur.cluster.i <- clusterAssign[i]
      Wi <- data[,,i]
      logdetWi <- logdetW[i]
      
      if(clusterSizes[cur.cluster.i] > 1){
        
        ## i is a non-singleton: |C_-i| = nClusters
        c.counts.noi <- clusterSizes
        c.counts.noi[cur.cluster.i] <- c.counts.noi[cur.cluster.i] - 1
        
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          ni <- as.numeric(c.counts.noi[x])
          if(ni <= 0) return(-Inf)
          
          ## Si = S_{x,-i}
          Si <- Sarr[,,x]
          if(x == cur.cluster.i){
            Si <- Si - Wi
          }
          
          log_prior <- log(GAMMA + ni)
          
          ## use the candidate existing cluster's own nu_x.
          log_like  <- log_pred_existing(Wi, Si, ni, nu[x],
                                         Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        ### m auxiliary new-cluster candidates
        aux_nu <- runif(m_aux, min = nu_L, max = nu_U)
        log_prior_new <- log(GAMMA) +
          (VN[nClusters+1] - VN[nClusters]) - log(m_aux)
        
        log_like_new <- sapply(1:m_aux, function(r){
          log_pred_new(Wi, aux_nu[r], Psi0, kappa0, logdetWi)
        })
        
        log_clusterProbs <- c(log_clusterProbs,
                              log_prior_new + log_like_new)
        
        ### sample among nClusters existing clusters and m_aux auxiliaries
        cluster.i <- sample(1:(nClusters+m_aux), size=1,
                            prob = exp(log_clusterProbs -
                                         sna::logSum(log_clusterProbs)))
        
        ### Update sufficient stats
        if(cluster.i == cur.cluster.i){
          #### stay: do nothing
          
        } else if(cluster.i <= nClusters){
          #### move to existing
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          clusterAssign[i] <- cluster.i
          
        } else {
          #### open new cluster
          #### retain the selected auxiliary nu as the new cluster's nu.
          aux.index <- cluster.i - nClusters
          
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          
          Snew <- array(0, dim=c(p,p,nClusters+1))
          Snew[,,1:nClusters] <- Sarr
          Snew[,,nClusters+1] <- Wi
          Sarr <- Snew
          
          A <- c(A, logdetWi)
          nu <- c(nu, aux_nu[aux.index]) 
          
          clusterAssign[i] <- nClusters + 1
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
        
      } else {
        
        ## i is a singleton: |C_-i| = nClusters - 1
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          if(x == cur.cluster.i) return(-Inf)
          
          ni <- as.numeric(clusterSizes[x])  # n_{x,-i}=n_x
          Si <- Sarr[,,x]
          
          log_prior <- log(GAMMA + ni)
          
          ## use the candidate existing cluster's own nu_x.
          log_like  <- log_pred_existing(Wi, Si, ni, nu[x],
                                         Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        ### Neal Algorithm 8 singleton update.
        ### The old singleton nu must be retained as one auxiliary candidate;
        ### only the remaining m_aux - 1 values are newly sampled from H_nu.
        aux_nu <- c(nu[cur.cluster.i],
                    runif(m_aux - 1, min = nu_L, max = nu_U))
        
        if(nClusters == 1){
          ## Only relevant when n = 1. There are no existing alternatives,
          ## so the common MFM new-cluster coefficient cancels.
          log_prior_new <- -log(m_aux)
        } else {
          log_prior_new <- log(GAMMA) +
            (VN[nClusters] - VN[nClusters-1]) - log(m_aux)
        }
        
        log_like_new <- sapply(1:m_aux, function(r){
          log_pred_new(Wi, aux_nu[r], Psi0, kappa0, logdetWi)
        })
        
        log_clusterProbs <- c(log_clusterProbs,
                              log_prior_new + log_like_new)
        
        cluster.i <- sample(1:(nClusters+m_aux), size=1,
                            prob = exp(log_clusterProbs -
                                         sna::logSum(log_clusterProbs)))
        
        if(cluster.i > nClusters){
          #### remain a singleton, but possibly replace the old nu
          #### by one of the fresh auxiliary values. If auxiliary 1 is chosen,
          #### both the partition and nu remain unchanged.
          aux.index <- cluster.i - nClusters
          nu[cur.cluster.i] <- aux_nu[aux.index]
          clusterAssign[i] <- cur.cluster.i
          
        } else {
          #### move i into existing cluster, then delete empty cur cluster and relabel
          #### update target stats first (indices are still old)
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          
          ## drop cur.cluster.i slice from Sarr/A
          Sarr <- Sarr[,,-cur.cluster.i, drop=FALSE]
          A    <- A[-cur.cluster.i]
          nu   <- nu[-cur.cluster.i] 
          
          ## relabel clusterAssign
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i,
                                  clusterAssign-1, clusterAssign)
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
      }
    } ## end i loop
    
    ## update each occupied cluster's nu_c separately by RW-MH,
    ## with Sigma_c analytically integrated out.
    nvec <- as.numeric(table(factor(clusterAssign,
                                    levels=1:max(clusterAssign))))
    nvec <- nvec[nvec > 0]
    nClusters <- length(nvec)
    
    for(c in 1:nClusters){
      nu_prop <- nu[c] + rnorm(1, mean=0, sd=proposal_sd)
      
      if(nu_prop >= nu_L && nu_prop <= nu_U){
        
        log_old <- log_target_nu_cluster(nu[c], Sarr[,,c], A[c], nvec[c],
                                         Psi0, kappa0)
        log_new <- log_target_nu_cluster(nu_prop, Sarr[,,c], A[c], nvec[c],
                                         Psi0, kappa0)
        
        log_alpha <- log_new - log_old
        if(log(runif(1)) < log_alpha){
          nu[c] <- nu_prop
        }
      }## else: reject
    }
    
    ## store the cluster-specific nu vector.
    History[[iter]] <- list(zout = clusterAssign, nuout = nu)
    
  } # end iter loop
  
  list(Iterates = History, VN = VN)
}



#########################################
DPM_Wishart_collapsed <- function(
    data, niterations,
    GAMMA,            # DP concentration (i.e., alpha)
    Psi0, kappa0,     # IW prior hyperparams
    nu_init, nu_L, nu_U,  # nu prior bounds + init
    initNClusters,   # number of initial clusters
    proposal_sd = 1, # RW-MH sd for each cluster-specific nu_c
    m_aux = 1        # number of Neal Algorithm 8 auxiliary nu values
){
  
  ## data: p x p x n array, each slice SPD
  n <- dim(data)[3]
  p <- nrow(data[,,1])
  
  if(m_aux < 1 || m_aux != as.integer(m_aux)){
    stop("m_aux must be a positive integer.")
  }
  m_aux <- as.integer(m_aux)
  
  ## Precompute logdet(Wi)(for updating z and nu)
  logdetW <- rep(NA_real_, n)
  for(i in 1:n){
    logdetW[i] <- logdet_spd(data[,,i])
  }
  
  ## initialize clusterAssign
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  ## for collapsed sampling
  ## For fast computation, we keep the sufficient stats in each cluster
  ##                S_c = sum W_i,  A_c = sum log|W_i|
  clusterSizes <- table(as.factor(clusterAssign))
  nClusters <- length(clusterSizes)
  
  Sarr <- array(0, dim=c(p,p,nClusters))
  A    <- rep(0, nClusters)
  for(i in 1:n){
    c <- clusterAssign[i]
    Sarr[,,c] <- Sarr[,,c] + data[,,i]
    A[c] <- A[c] + logdetW[i]
  }
  
  ## A scalar nu_init is replicated; alternatively pass a vector of length initNClusters.
  if(length(nu_init) == 1){
    nu <- rep(nu_init, nClusters)
  } else {
    if(length(nu_init) != nClusters){
      stop("nu_init must be a scalar or a vector of length initNClusters.")
    }
    nu <- nu_init
  }
  
  History <- vector("list", niterations)
  
  ## start MCMC
  for(iter in 1:niterations){
    
    #if(iter %% 1000 == 0){
    #  print(iter)
    #}
    
    ## update z
    clusterSizes <- table(as.factor(clusterAssign))
    nClusters <- length(clusterSizes)
    
    for(i in 1:n){
      
      cur.cluster.i <- clusterAssign[i]
      Wi <- data[,,i]
      logdetWi <- logdetW[i]
      
      if(clusterSizes[cur.cluster.i] > 1){
        
        ## i is a non-singleton: |C_-i| = nClusters
        c.counts.noi <- clusterSizes
        c.counts.noi[cur.cluster.i] <- c.counts.noi[cur.cluster.i] - 1
        
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          ni <- as.numeric(c.counts.noi[x])
          if(ni <= 0) return(-Inf)
          
          ### Si = S_{x,-i}
          Si <- Sarr[,,x]
          if(x == cur.cluster.i){
            Si <- Si - Wi
          }
          
          log_prior <- log(ni)
          
          ## use the candidate existing cluster's own nu_x.
          log_like  <- log_pred_existing(Wi, Si, ni, nu[x],
                                         Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        ### m auxiliary new-cluster candidates
        aux_nu <- runif(m_aux, min = nu_L, max = nu_U)
        log_prior_new <- log(GAMMA) - log(m_aux)
        
        log_like_new <- sapply(1:m_aux, function(r){
          log_pred_new(Wi, aux_nu[r], Psi0, kappa0, logdetWi)
        })
        
        log_clusterProbs <- c(log_clusterProbs,
                              log_prior_new + log_like_new)
        
        ## sample among nClusters existing clusters and m_aux auxiliaries
        cluster.i <- sample(1:(nClusters+m_aux), size=1,
                            prob = exp(log_clusterProbs -
                                         sna::logSum(log_clusterProbs)))
        
        ## update sufficient stats
        if(cluster.i == cur.cluster.i){
          ### stay: do nothing
          
        } else if(cluster.i <= nClusters){
          ### move to existing
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          clusterAssign[i] <- cluster.i
          
        } else {
          ### open new cluster
          ### retain the selected auxiliary nu as the new cluster's nu.
          aux.index <- cluster.i - nClusters
          
          Sarr[,,cur.cluster.i] <- Sarr[,,cur.cluster.i] - Wi
          A[cur.cluster.i] <- A[cur.cluster.i] - logdetWi
          
          Snew <- array(0, dim=c(p,p,nClusters+1))
          Snew[,,1:nClusters] <- Sarr
          Snew[,,nClusters+1] <- Wi
          Sarr <- Snew
          
          A <- c(A, logdetWi)
          nu <- c(nu, aux_nu[aux.index])
          
          clusterAssign[i] <- nClusters + 1
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
        
      } else {
        
        ## i is a singleton: |C_-i| = nClusters - 1
        ### existing clusters
        log_clusterProbs <- sapply(1:nClusters, function(x){
          
          if(x == cur.cluster.i) return(-Inf)
          
          ni <- as.numeric(clusterSizes[x])  # n_{x,-i}=n_x
          Si <- Sarr[,,x]
          
          log_prior <- log(ni)
          
          ## use the candidate existing cluster's own nu_x.
          log_like  <- log_pred_existing(Wi, Si, ni, nu[x],
                                         Psi0, kappa0, logdetWi)
          log_prior + log_like
        })
        
        
        ### Neal Algorithm 8 singleton update.
        ### The old singleton nu must be retained as one auxiliary candidate;
        ### only the remaining m_aux - 1 values are newly sampled from H_nu.
        aux_nu <- c(nu[cur.cluster.i],
                    runif(m_aux - 1, min = nu_L, max = nu_U))
        
        log_prior_new <- log(GAMMA) - log(m_aux)
        
        log_like_new <- sapply(1:m_aux, function(r){
          log_pred_new(Wi, aux_nu[r], Psi0, kappa0, logdetWi)
        })
        
        log_clusterProbs <- c(log_clusterProbs,
                              log_prior_new + log_like_new)
        
        cluster.i <- sample(1:(nClusters+m_aux), size=1,
                            prob = exp(log_clusterProbs -
                                         sna::logSum(log_clusterProbs)))
        
        if(cluster.i > nClusters){
          ## remain a singleton, but possibly replace the old nu
          ## by one of the fresh auxiliary values. If auxiliary 1 is chosen,
          ## both the partition and nu remain unchanged.
          aux.index <- cluster.i - nClusters
          nu[cur.cluster.i] <- aux_nu[aux.index]
          clusterAssign[i] <- cur.cluster.i
          
        } else {
          ## move i into existing cluster, then delete empty cur cluster and relabel
          ## update target stats first (indices are still old)
          Sarr[,,cluster.i] <- Sarr[,,cluster.i] + Wi
          A[cluster.i] <- A[cluster.i] + logdetWi
          
          ## drop cur.cluster.i slice from Sarr/A
          Sarr <- Sarr[,,-cur.cluster.i, drop=FALSE]
          A    <- A[-cur.cluster.i]
          nu   <- nu[-cur.cluster.i]  
          
          ## relabel clusterAssign
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i,
                                  clusterAssign-1, clusterAssign)
        }
        
        clusterSizes <- table(as.factor(clusterAssign))
        nClusters <- length(clusterSizes)
      }
    } # end i loop
    
    ## update each occupied cluster's nu_c separately by RW-MH,
    ## with Sigma_c analytically integrated out.
    nvec <- as.numeric(table(factor(clusterAssign,
                                    levels=1:max(clusterAssign))))
    nvec <- nvec[nvec > 0]
    nClusters <- length(nvec)
    
    for(c in 1:nClusters){
      nu_prop <- nu[c] + rnorm(1, mean=0, sd=proposal_sd)
      
      if(nu_prop >= nu_L && nu_prop <= nu_U){
        
        log_old <- log_target_nu_cluster(nu[c], Sarr[,,c], A[c], nvec[c],
                                         Psi0, kappa0)
        log_new <- log_target_nu_cluster(nu_prop, Sarr[,,c], A[c], nvec[c],
                                         Psi0, kappa0)
        
        log_alpha <- log_new - log_old
        if(log(runif(1)) < log_alpha){
          nu[c] <- nu_prop
        }
      }## else: reject
    }
    
    ## store the cluster-specific nu vector.
    History[[iter]] <- list(zout = clusterAssign, nuout = nu)
    
  } # end iter loop
  
  list(Iterates = History)
}

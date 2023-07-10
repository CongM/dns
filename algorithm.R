#######################################################
##                                                   ##
## Supplementary Code                                ##
##                                                   ##
## Cong Mu, Youngser Park, and Carey E. Priebe       ##
##                                                   ##
#######################################################




#### Required packages
library(ggplot2)
library(dplyr)
library(irlba)
library(mclust)
library(igraph)
library(grdpg)

source("codeFunc_approx_chernoff_optimize_v2.R")


#### Some utility functions
ASE_GMM <- function(A, dmax = 10, dhat = NULL, G = 3:10, seed = 2020) {
  diag(A) <- rowSums(A) / (nrow(A)-1)
  embedding <- irlba(A, dmax)
  s <- embedding$d
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
  Xhat <- embedding$u[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  
  set.seed(seed)
  model <- Mclust(Xhat, G)
  tauhat <- model$classification
  
  return(tauhat)
}

AtoE <- function(A) {
  A[lower.tri(A)] <- 0
  E <- which(A==1, arr.ind = TRUE)
  return(E)
}

EtoA <- function(n, E) {
  A <- sparseMatrix(i = E[,1], j = E[,2], x = 1, dims = c(n,n), symmetric = TRUE)
  return(A)
}

Chernoffactive <- function(ij, xihat, kstar) {
  iblock <- xihat[ij[1]]
  jblock <- xihat[ij[2]]
  indi <- iblock %in% kstar
  indj <- jblock %in% kstar
  if (indi & indj) {
    return(1)
  } else {
    return(0)
  }
}


#### Algorithm 1
Algo1 <- function(n, E, ind0, p1, dmax = 10, dhat = NULL, G = 3:10, seed = 2020) {
  ## Step 1
  indtemp <- setdiff(1:nrow(E), ind0)
  ind1 <- sample(indtemp, floor(nrow(E)*p1))
  
  ## Step 2
  dE <- E[c(ind0,ind1),]
  A <- EtoA(n, dE)
  
  ## Step 3
  diag(A) <- rowSums(A) / (nrow(A)-1)
  embedding <- irlba(A, dmax)
  s <- embedding$d
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
  Xhat <- embedding$u[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  
  ## Step 4
  set.seed(seed)
  model <- Mclust(Xhat, G)
  tauhat <- model$classification
  
  return(tauhat)
}


#### Algorithm 2 
Algo2 <- function(n, E, ind0, p1, dmax = 10, dhat = NULL, G = 3:10, seed = 2020) {
  ## Step 1
  dE <- E[ind0,]
  A <- EtoA(n, dE)
  
  ## Step 2
  diag(A) <- rowSums(A) / (nrow(A)-1)
  embedding <- irlba(A, dmax)
  s <- embedding$d
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
  Xhat <- embedding$u[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  
  ## Step 3
  set.seed(seed)
  model <- Mclust(Xhat, G)
  xihat <- model$classification
  
  ## Step 4
  Khat <- model$G
  pihat <- rep(0, Khat)
  for (k in 1:Khat) {
    pihat[k] <- sum(xihat==k) / n
  }
  
  ## Step 5
  Ipq <- getIpq(A, dhat)
  muhats <- matrix(model$parameters$mean, nrow = dhat)
  Bhat <- t(muhats) %*% Ipq %*% muhats
  
  ## Step 6
  CA <- approx_chernoff_opt_v2(Bhat, pihat)
  kstar <- CA[2:3]
  
  ## Step 7
  indtemp <- setdiff(1:nrow(E), ind0)
  Etemp <- E[indtemp,]
  tempind <- apply(Etemp, 1, Chernoffactive, xihat, kstar)
  indstar <- which(tempind==1)
  indtemp1 <- indtemp[indstar]
  n1 <- floor(floor(nrow(E)*p1)*sum(pihat[kstar])^2)
  if (n1 < length(indtemp1)) {
    ind1 <- sample(indtemp1, n1)
  } else {
    ind1 <- indtemp1
  }
  n2 <- floor(nrow(E)*p1) - length(ind1)
  indtemp2 <- setdiff(indtemp, ind1)
  ind2 <- sample(indtemp2, n2)
  
  ## Step 8
  dE <- E[c(ind0,ind1,ind2),]
  A <- EtoA(n, dE)
  
  ## Step 9
  diag(A) <- rowSums(A) / (nrow(A)-1)
  embedding <- irlba(A, dmax)
  s <- embedding$d
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
  Xhat <- embedding$u[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  
  ## Step 10
  set.seed(seed)
  model <- Mclust(Xhat, G)
  tauhat <- model$classification
  
  return(tauhat)
}




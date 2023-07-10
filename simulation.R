#######################################################
##                                                   ##
## Supplementary Code                                ##
##                                                   ##
## Cong Mu, Youngser Park, and Carey E. Priebe       ##
##                                                   ##
#######################################################




source("algorithm.R")


#### Figure 4
K <- 4

## Balanced case
pi <- rep(1/K, K)

# ## Unbalanced case
# pi <- c(1/8, 1/8, 3/8, 3/8)

q1 <- 0.2
q2 <- 0.4
q3 <- 0.5
q4 <- 0.9
nu <- c(q1, q2, q3, q4)
B <- nu %*% t(nu)

p0 <- 0.01
B0 <- p0*B

result0 <- approx_chernoff_opt_v2(B0, pi)
rhoB0 <- result0[1]
kstar <- result0[2:3]

onestar <- matrix(0, nrow = nrow(B), ncol = ncol(B))
onestar[kstar,kstar] <- 1

## Balanced case
p1star <- 0.1940

# ## Unbalanced case
# p1star <- 0.2201

B1tilde <- B0 + (p1star/(sum(pi[kstar])^2)) * (B*onestar)
result1tilde <- approx_chernoff_opt_v2(B1tilde, pi)
kstar1 <- result1tilde[2:3]

p11max <- 1 - p0 - (p1star/(sum(pi[kstar])^2)) + p1star

dmax <- 10
dhat <- 1
G <- K
seed <- 2020

n <- 12000
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

Aseeds <- 1:50

## Balanced case
p1s <- c(seq(from = 0, to = p1star, length.out = 10), seq(from = p1star, to = p11max, length.out = 11)[-1])

# ## Unbalanced case
# p1s <- c(seq(from = 0, to = p1star, length.out = 12), seq(from = p1star, to = p11max, length.out = 9)[-1])

allresults <- data.frame()
for (p1 in p1s) {
  cat(p1, "\n")
  B1 <- B0 + p1*B
  if (p1 < p1star) {
    B1tilde <- B0 + (p1/(sum(pi[kstar])^2)) * (B*onestar)
  } else {
    B1tilde <- B0 + (p1-p1star)*B + (p1star/(sum(pi[kstar])^2)) * (B*onestar)
  }
  
  ARI1 <- c()
  ARI1tilde <- c()
  elapsed1 <- c()
  elapsed1tilde <- c()
  for (Aseed in Aseeds) {
    # cat(Aseed, "\n")
    
    set.seed(Aseed)
    g1 <- sample_sbm(n, pref.matrix = B1, block.sizes = block_size)
    A1 <- g1[]
    
    set.seed(Aseed)
    g1tilde <- sample_sbm(n, pref.matrix = B1tilde, block.sizes = block_size)
    A1tilde <- g1tilde[]
    
    start_time1 <- Sys.time()
    tauhat1 <- ASE_GMM(A1, dmax, dhat, G, seed)
    end_time1 <- Sys.time()
    elapsed1 <- c(elapsed1, as.numeric(difftime(end_time1,start_time1,units="secs")))
    ARI1 <- c(ARI1, adjustedRandIndex(blocks, tauhat1))
    
    start_time1tilde <- Sys.time()
    tauhat1tilde <- ASE_GMM(A1tilde, dmax, dhat, G, seed)
    end_time1tilde <- Sys.time()
    elapsed1tilde <- c(elapsed1tilde, as.numeric(difftime(end_time1tilde,start_time1tilde,units="secs")))
    ARI1tilde <- c(ARI1tilde, adjustedRandIndex(blocks, tauhat1tilde))
  }
  
  tempresults1 <- data.frame(p1, ARI = ARI1, Elapsed = elapsed1, Scheme = "B1", Aseed = Aseeds)
  tempresults2 <- data.frame(p1, ARI = ARI1tilde, Elapsed = elapsed1tilde, Scheme = "B1tilde", Aseed = Aseeds)
  allresults <- rbind(allresults, tempresults1, tempresults2)
}

dat <- allresults %>%
  group_by(p1, Scheme) %>%
  summarise(ARI_mean = mean(ARI), 
            ARI_se = sd(ARI) / sqrt(max(Aseeds)),
            Elapsed_mean = mean(Elapsed), 
            Elapsed_se = sd(Elapsed) / sqrt(max(Aseeds)))

col <- c("#F8766D", "#00BFC4")
ggplot(dat) +
  theme_bw() +
  geom_line(aes(x=p1,y=ARI_mean,color=Scheme), size = 2) +
  geom_ribbon(aes(x=p1,ymin=ARI_mean-ARI_se,ymax=ARI_mean+ARI_se,fill=Scheme), alpha = 0.3) +
  geom_hline(yintercept = 1, color = "black", size = 3) +
  geom_vline(xintercept = p1star, color = "black", linetype = 2, size = 1) +
  labs(x = expression(p[1]), y = "ARI", color = " ", fill = " ") +
  scale_color_manual(labels = c(expression(B[1]), expression(paste(widetilde(B)[1]," / ",widetilde(B)[1],"*"))), values = col) +
  scale_fill_manual(labels = c(expression(B[1]), expression(paste(widetilde(B)[1]," / ",widetilde(B)[1],"*"))), values = col) +
  theme(axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 30),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 30),
        legend.title = element_text(size = 32),
        legend.text = element_text(size = 30),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")


#### Figure 5
addCovariates <- FALSE
dmax <- 10
G <- 3:10
seed <- 2020

n <- 4000
K <- 4
d <- 1

q1 <- 0.2
q2 <- 0.4
q3 <- 0.5
q4 <- 0.9
latent <- cbind(q1, q2, q3, q4)

## Balanced case
pi <- rep(1/K, K)

# ## Unbalanced case
# pi <- c(1/8, 1/8, 3/8, 3/8)

block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)

p0 <- 0.15
Aseeds <- 1:50
Es <- list()
ind0s <- list()
for (Aseed in Aseeds) {
  tempA <- generateA(n, P, seed = Aseed)
  diag(tempA) <- 0
  tempE <- AtoE(tempA)
  tempind0 <- sample(nrow(tempE), floor(nrow(tempE)*p0))
  Es[[Aseed]] <- tempE
  ind0s[[Aseed]] <- tempind0
}

dhat <- rankMatrix(B)

p1s <- seq(from = 0, to = 1-p0, length.out = 20)
allresults_cg <- data.frame()
allresults_gt <- data.frame()
for (p1 in p1s) {
  cat(p1, "\n")
  ARI1_cg <- c()
  ARI2_cg <- c()
  ARI1_gt <- c()
  ARI2_gt <- c()
  ARI12 <- c()
  elapsed1 <- c()
  elapsed2 <- c()
  for (Aseed in Aseeds) {
    # cat(Aseed, "\n")
    
    E <- Es[[Aseed]]
    ind0 <- ind0s[[Aseed]]
    
    start_time1 <- Sys.time()
    tauhat1 <- Algo1(n, E, ind0, p1, dmax, dhat, G, seed)
    end_time1 <- Sys.time()
    elapsed1 <- c(elapsed1, as.numeric(difftime(end_time1,start_time1,units="secs")))
    
    start_time2 <- Sys.time()
    tauhat2 <- Algo2(n, E, ind0, p1, dmax, dhat, G, seed)
    end_time2 <- Sys.time()
    elapsed2 <- c(elapsed2, as.numeric(difftime(end_time2,start_time2,units="secs")))
    
    A <- EtoA(n, E)
    diag(A) <- rowSums(A) / (nrow(A)-1)
    embedding <- irlba(A, dmax)
    s <- embedding$d
    Xhat <- embedding$u[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
    set.seed(seed)
    model <- Mclust(Xhat, G)
    tauhat <- model$classification
    
    ARI1_cg <- c(ARI1_cg, adjustedRandIndex(tauhat, tauhat1))
    ARI2_cg <- c(ARI2_cg, adjustedRandIndex(tauhat, tauhat2))
    ARI12 <- c(ARI12, adjustedRandIndex(tauhat1, tauhat2))
    ARI1_gt <- c(ARI1_gt, adjustedRandIndex(blocks, tauhat1))
    ARI2_gt <- c(ARI2_gt, adjustedRandIndex(blocks, tauhat2))
  }
  tempresults1_cg <- data.frame(p1, ARI = ARI1_cg, Algorithm = "Algorithm 1", Aseed = Aseeds, elapsed = elapsed1)
  tempresults2_cg <- data.frame(p1, ARI = ARI2_cg, Algorithm = "Algorithm 2", Aseed = Aseeds, elapsed = elapsed2)
  tempresults12 <- data.frame(p1, ARI = ARI12, Algorithm = "Algorithm 12", Aseed = Aseeds, elapsed = 0)
  allresults_cg <- rbind(allresults_cg, tempresults1_cg, tempresults2_cg, tempresults12)
  tempresults1_gt <- data.frame(p1, ARI = ARI1_gt, Algorithm = "Algorithm 1", Aseed = Aseeds, elapsed = elapsed1)
  tempresults2_gt <- data.frame(p1, ARI = ARI2_gt, Algorithm = "Algorithm 2", Aseed = Aseeds, elapsed = elapsed2)
  allresults_gt <- rbind(allresults_gt, tempresults1_gt, tempresults2_gt)
}

dat <- allresults_gt %>%
  group_by(p1, Algorithm) %>%
  summarise(ARI_mean = mean(ARI), 
            ARI_se = sd(ARI) / sqrt(max(Aseeds)),
            elapsed_mean = mean(elapsed),
            elapsed_se = sd(elapsed) / sqrt(max(Aseeds)))

ggplot(dat) +
  theme_bw() +
  geom_line(aes(x=p1,y=ARI_mean,color=Algorithm), size = 2) +
  geom_ribbon(aes(x=p1,ymin=ARI_mean-ARI_se,ymax=ARI_mean+ARI_se,fill=Algorithm), alpha = 0.3) +
  geom_hline(yintercept = 1, color = "black", size = 3) +
  labs(x = expression(p[1]), y = "ARI", color = " ", fill = " ") +
  theme(axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 30),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 30),
        legend.title = element_text(size = 32),
        legend.text = element_text(size = 30),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")




#######################################################
##                                                   ##
## Supplementary Code                                ##
##                                                   ##
## Cong Mu, Youngser Park, and Carey E. Priebe       ##
##                                                   ##
#######################################################




source("algorithm.R")


#### Figure 8
SNDresults <- data.frame()

# ## LastFM
# Dataset <- "LastFM"
# el <- read.csv("./lasftm_asia/lastfm_asia_edges.csv") %>% 
#   mutate(node_1 = node_1 + 1, node_2 = node_2 + 1) %>%
#   as.matrix()
# g <- graph_from_edgelist(el, directed = FALSE)
# A <- g[]
# diag(A) <- 0
# 
# blocks <- read.csv("./lasftm_asia/lastfm_asia_target.csv") %>%
#   mutate(id = id + 1, target = target + 1)
# tau <- blocks$target
# 
# p0 <- 0.15
# G <- 3:20

## Facebook
Dataset <- "Facebook"
el <- read.csv("./facebook_large/musae_facebook_edges.csv") %>% 
  mutate(id_1 = id_1 + 1, id_2 = id_2 + 1) %>%
  as.matrix()
g <- graph_from_edgelist(el, directed = FALSE)
A <- g[]
diag(A) <- 0

blocks <- read.csv("./facebook_large/musae_facebook_target.csv") %>%
  mutate(id = id + 1,
         block = ifelse(page_type == "company", 1, 
                        ifelse(page_type == "government", 2, 
                               ifelse(page_type == "politician", 3, 4))))
tau <- blocks$block

p0 <- 0.35
G <- 3:10


n <- nrow(A)
E <- AtoE(A)

dhat <- NULL
dmax <- 30
seeds <- 1:100

set.seed(2020)
ind0 <- sample(nrow(E), floor(nrow(E)*p0))

diag(A) <- rowSums(A) / (nrow(A)-1)
embedding <- irlba(A, dmax)
s <- embedding$d
dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
Xhat <- embedding$u[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
set.seed(2020)
model <- Mclust(Xhat, G)
tauhat <- model$classification

p1s <- seq(from = 0.05, length.out = 5, by = 0.05)
for (p1 in p1s) {
  cat(p1, "\n")
  
  ARI_Algo1 <- c()
  ARI_Algo2 <- c()
  ARI_Algo12 <- c()
  ARI_Algo1_cg <- c()
  ARI_Algo2_cg <- c()
  elapsed1 <- c()
  elapsed2 <- c()
  for (seed in seeds) {
    # cat(seed, "\n")
    set.seed(seed)
    
    start_time1 <- Sys.time()
    tauhat1 <- Algo1(n, E, ind0, p1, dmax, dhat, G, seed)
    end_time1 <- Sys.time()
    elapsed1 <- c(elapsed1, as.numeric(difftime(end_time1,start_time1,units="secs")))
    
    start_time2 <- Sys.time()
    tauhat2 <- Algo2(n, E, ind0, p1, dmax, dhat, G, seed)
    end_time2 <- Sys.time()
    elapsed2 <- c(elapsed2, as.numeric(difftime(end_time2,start_time2,units="secs")))
    
    ARI_Algo1 <- c(ARI_Algo1, adjustedRandIndex(tau, tauhat1))
    ARI_Algo2 <- c(ARI_Algo2, adjustedRandIndex(tau, tauhat2))
    
    ARI_Algo1_cg <- c(ARI_Algo1_cg, adjustedRandIndex(tauhat, tauhat1))
    ARI_Algo2_cg <- c(ARI_Algo2_cg, adjustedRandIndex(tauhat, tauhat2))
    ARI_Algo12 <- c(ARI_Algo12, adjustedRandIndex(tauhat1, tauhat2))
  }
  
  tempresults1 <- data.frame(p0, p1, ARI = ARI_Algo1, Algorithm = "Algorithm 1", Dataset, Block = "GT", elapsed = elapsed1, Seed = seeds)
  tempresults2 <- data.frame(p0, p1, ARI = ARI_Algo2, Algorithm = "Algorithm 2", Dataset, Block = "GT", elapsed = elapsed2, Seed = seeds)
  tempresults1_cg <- data.frame(p0, p1, ARI = ARI_Algo1_cg, Algorithm = "Algorithm 1", Dataset, Block = "CG", elapsed = elapsed1, Seed = seeds)
  tempresults2_cg <- data.frame(p0, p1, ARI = ARI_Algo2_cg, Algorithm = "Algorithm 2", Dataset, Block = "CG", elapsed = elapsed2, Seed = seeds)
  tempresults12 <- data.frame(p0, p1, ARI = ARI_Algo12, Algorithm = "Algorithm 12", Dataset, Block = "Algo12", elapsed = 0, Seed = seeds)
  SNDresults <- rbind(SNDresults, tempresults1, tempresults2, tempresults1_cg, tempresults2_cg, tempresults12)
}

SNDresults <- SNDresults %>%
  arrange(Dataset, Block, p0, p1, Algorithm, Seed)

# ## LastFM
# dat <- SNDresults %>%
#   filter(Algorithm != "Algorithm 12" & Block == "GT" & p0 == 0.15) %>%
#   group_by(Dataset, p0, p1, Algorithm) %>%
#   summarise(ARI_mean = mean(ARI), 
#             ARI_se = sd(ARI) / sqrt(max(seeds)),
#             elapsed_mean = mean(elapsed),
#             elapsed_se = sd(elapsed) / sqrt(max(seeds)))

## Facebook
dat <- SNDresults %>%
  filter(Algorithm != "Algorithm 12" & Block == "GT" & p0 == 0.35) %>%
  group_by(Dataset, p0, p1, Algorithm) %>%
  summarise(ARI_mean = mean(ARI), 
            ARI_se = sd(ARI) / sqrt(max(seeds)),
            elapsed_mean = mean(elapsed),
            elapsed_se = sd(elapsed) / sqrt(max(seeds)))

ggplot(dat) +
  theme_bw() +
  geom_line(aes(x=p1,y=ARI_mean,color=Algorithm), size = 2) +
  geom_ribbon(aes(x=p1,ymin=ARI_mean-ARI_se,ymax=ARI_mean+ARI_se,fill=Algorithm), alpha = 0.3) +
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




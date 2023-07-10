#######################################################
##                                                   ##
## Supplementary Code                                ##
##                                                   ##
## Cong Mu, Youngser Park, and Carey E. Priebe       ##
##                                                   ##
#######################################################




source("algorithm.R")


#### Some parameters
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

result <- approx_chernoff_opt_v2(B, pi)
rhoB <- result[1]


#### Figure 1
p0 <- 0.01
B0 <- p0*B

result0 <- approx_chernoff_opt_v2(B0, pi)
rhoB0 <- result0[1]
kstar <- result0[2:3]

onestar <- matrix(0, nrow = nrow(B), ncol = ncol(B))
onestar[kstar,kstar] <- 1

p0s <- seq(from = 0.0001, to = 0.9999, by = 0.0001)
rhoB0s <- c()
for (p0 in p0s) {
  B0 <- p0*B
  result0 <- approx_chernoff_opt_v2(B0, pi)
  rhoB0s <- c(rhoB0s, result0[1])
}

allresults <- rbind(data.frame(p = p0s, rho = rhoB, scheme = "B"),
                    data.frame(p = p0s, rho = rhoB0s, scheme = "pB"))

ggplot(allresults) +
  theme_bw() +
  geom_line(aes(x=p,y=rho,color=scheme), size = 2) +
  labs(y = expression(rho), color = " ") +
  theme(axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 30),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 30),
        legend.title = element_text(size = 32),
        legend.text = element_text(size = 30),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")


#### Figure 2
p0 <- 0.01
B0 <- p0*B

result0 <- approx_chernoff_opt_v2(B0, pi)
rhoB0 <- result0[1]
kstar <- result0[2:3]

onestar <- matrix(0, nrow = nrow(B), ncol = ncol(B))
onestar[kstar,kstar] <- 1

p1max <- (1-p0) * sum(pi[kstar])^2

p1s <- seq(from = 0, to = 1-p0, by = 0.0001)
rhoB1s <- c()
rhoB1tildes <- c()
klstars <- data.frame()
for (p1 in p1s) {
  B1 <- B0 + p1*B
  result1 <- approx_chernoff_opt_v2(B1, pi)
  rhoB1s <- c(rhoB1s, result1[1])
  
  if (p1 <= p1max) {
    B1tilde <- B0 + (p1/(sum(pi[kstar])^2)) * (B*onestar)
    result1tilde <- approx_chernoff_opt_v2(B1tilde, pi)
    rhoB1tildes <- c(rhoB1tildes, result1tilde[1])
    klstars <- rbind(klstars, data.frame(p1, result1tilde[2], result1tilde[3]))
  }
}

row.names(klstars) <- NULL
names(klstars)[2:3] <- c("kstar", "lstar")

## Balanced case
indstar <- 1940

# ## Unbalanced case
# indstar <- 2201

allresults <- rbind(data.frame(p1 = p1s[1:indstar], rho = rhoB, scheme = "B"),
                    data.frame(p1 = p1s[1:indstar], rho = rhoB0, scheme = "B0"),
                    data.frame(p1 = p1s[1:indstar], rho = rhoB1s[1:indstar], scheme = "B1"),
                    data.frame(p1 = p1s[1:indstar], rho = rhoB1tildes[1:indstar], scheme = "B1tilde"))

col <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")

## Balanced case
ggplot(allresults) +
  theme_bw() +
  xlim(0, 0.20) + 
  geom_line(aes(x=p1,y=rho,color=scheme), size = 2) +
  labs(x = expression(p[1]), y = expression(rho), color = " ") +
  scale_color_manual(labels = c('B', expression(B[0]), expression(B[1]), expression(widetilde(B)[1])), values = col) +
  theme(axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 30),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 30),
        legend.title = element_text(size = 32),
        legend.text = element_text(size = 30),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")

# ## Unbalanced case
# ggplot(allresults) +
#   theme_bw() +
#   geom_line(aes(x=p1,y=rho,color=scheme), size = 2) +
#   labs(x = expression(p[1]), y = expression(rho), color = " ") +
#   scale_color_manual(labels = c('B', expression(B[0]), expression(B[1]), expression(widetilde(B)[1])), values = col) +
#   theme(axis.title.x = element_text(size = 32),
#         axis.text.x = element_text(size = 30),
#         axis.title.y = element_text(size = 32),
#         axis.text.y = element_text(size = 30),
#         legend.title = element_text(size = 32),
#         legend.text = element_text(size = 30),
#         legend.key.width = unit(2.5, "cm"),
#         legend.key.height = unit(1.5, "cm"),
#         legend.position = "top")


#### Figure 3
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

p1s <- seq(from = p1star, to = p11max, by = 0.0001)
rhoB1s <- c()
rhoB1tildes <- c()
for (p1 in p1s) {
  B1 <- B0 + p1*B
  result1 <- approx_chernoff_opt_v2(B1, pi)
  rhoB1s <- c(rhoB1s, result1[1])
  
  B1tilde <- B0 + (p1-p1star)*B + (p1star/(sum(pi[kstar])^2)) * (B*onestar)
  result1tilde <- approx_chernoff_opt_v2(B1tilde, pi)
  rhoB1tildes <- c(rhoB1tildes, result1tilde[1])
}

allresults <- rbind(data.frame(p1 = p1s, rho = rhoB, scheme = "B"),
                    data.frame(p1 = p1s, rho = rhoB0, scheme = "B0"),
                    data.frame(p1 = p1s, rho = rhoB1s, scheme = "B1"),
                    data.frame(p1 = p1s, rho = rhoB1tildes, scheme = "B1tilde"))

col <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
ggplot(allresults) +
  theme_bw() +
  geom_line(aes(x=p1,y=rho,color=scheme), size = 2) +
  labs(x = expression(p[1]), y = expression(rho), color = " ") +
  scale_color_manual(labels = c('B', expression(B[0]), expression(B[1]), expression(paste(widetilde(B)[1],"*"))), values = col) +
  theme(axis.title.x = element_text(size = 32),
        axis.text.x = element_text(size = 30),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 30),
        legend.title = element_text(size = 32),
        legend.text = element_text(size = 30),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")





# Load packages and data --------------------------------------------------

library(tidyverse)
library(MCMCvis)
library(doSNOW)
library(foreach)
library(rjags)
library(MCMCvis)
library(ggthemes)

# Install with devtools::install_github("thomasp85/patchwork")
library(patchwork)

results_maps_hmsc <- readRDS("results_maps_hmsc.rds")
results_cjspop <- readRDS("results_cjspop.rds")
results_dem <- readRDS("results_dem.rds")
results_pR <- readRDS("results_pR.rds")
results_auc <- readRDS("results_auc.rds")

N <- results_dem$N
pR <- list()
for (i in 1:length(results_pR)) {
  
  pR[[i]] <- results_pR[[i]]
  if (any(pR[[i]]==1)) pR[[i]][which(pR[[i]]==1)] <- 0.9999
  
}

spcode <- results_cjspop$spcode
names(results_auc) <- spcode

R <- results_dem$R
popR <- list()
for (i in 1:17) {
  
  tempR <- c()
  for (h in 1:nrow(R[[i]]$R)) {
    
    index <- which(results_cjspop$chdata[[i]]$effort_year[h,]>0)
    y <- R[[i]]$R[h,index]
    tempR[h] <- prod(y)^(1/length(y))
    
  }
  
  popR[[i]] <- tempR
  
}
names(popR) <- spcode

# Load jags functions
source("models_jags.R")

# Can prob of occ predict dem params? -------------------------------------

x_std <- map(results_maps_hmsc, function(x) (x - mean(x))/sd(x))

## R vs Pr
RvsPr <- foreach (i=1:length(popR)) %do% {
  
  inits <- list(
    list(alpha = 20, beta = -20, sigma = 1),
    list(alpha = 0, beta = 0, sigma = 50),
    list(alpha = -20, beta = 20, sigma = 100)
  )
  
  res <- lm.jags(data = list(x = x_std[[i]], 
                             y = log(popR[[i]]), 
                             n = length(popR[[i]])),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 10,
                 seed_no = 19)
  
  return(res)
  
}

## N vs Pr
NvsPr <- foreach (i=1:length(popR)) %do% {
  
  inits <- list(
    list(alpha = 20, beta = -20, sigma = 1),
    list(alpha = 0, beta = 0, sigma = 50),
    list(alpha = -20, beta = 20, sigma = 100)
  )
  
  res <- lm.jags(data = list(x = x_std[[i]],
                             y = log(N[[i]]$popN), 
                             n = length(x_std[[i]])),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 10,
                 seed_no = 19)
  
  return(res)
  
}

## High N vs Pr
hNvsPr <- foreach (i=1:length(avgR)) %do% {
  
  inits <- list(
    list(alpha = 20, beta = -20, sigma = 1),
    list(alpha = 0, beta = 0, sigma = 50),
    list(alpha = -20, beta = 20, sigma = 100)
  )
  
  res <- lm.jags(data = list(x = x_std[[i]],
                             y = log(N[[i]]$highN), 
                             n = length(x_std[[i]])),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 10,
                 seed_no = 19)
  
  return(res)
  
}


## PR vs Pr
pRvsPr <- foreach (i=1:(length(popR))) %do% {
  
  inits <- list(
    list(alpha = 5, beta = -5),
    list(alpha = 0, beta = 0),
    list(alpha = -5, beta = 5)
  )
  
  res <- beta.lm.jags(data = list(x = x_std[[i]], y = pR[[i]], 
                                  n = length(x_std[[i]])),
                      inits = inits,
                      n.chains = 3,
                      n.adapt = 1000,
                      n.update = 5000,
                      n.iter = 5000,
                      n.thin = 1,
                      seed_no = 19)
  
  return(res)
  
}

## Check convergence
map_dbl(RvsPr, function(x) max(x$mcmc_sum[,"Rhat"]))
map_dbl(NvsPr, function(x) max(x$mcmc_sum[,"Rhat"]))
map_dbl(hNvsPr, function(x) max(x$mcmc_sum[,"Rhat"]))
map_dbl(pRvsPr, function(x) max(x$mcmc_sum[,"Rhat"]))

map_dbl(RvsPr, function(x) min(x$mcmc_sum[,"n.eff"]))
map_dbl(NvsPr, function(x) min(x$mcmc_sum[,"n.eff"]))
map_dbl(hNvsPr, function(x) min(x$mcmc_sum[,"n.eff"]))
map_dbl(pRvsPr, function(x) min(x$mcmc_sum[,"n.eff"]))


## Plot

beta_chains <-  list()
beta_chains[[1]] <- map(RvsPr, function(x) MCMCchains(x$mcmc_samples, params = "beta"))
beta_chains[[2]] <- map(NvsPr, function(x) MCMCchains(x$mcmc_samples, params = "beta"))
beta_chains[[3]] <- map(hNvsPr, function(x) MCMCchains(x$mcmc_samples, params = "beta"))
beta_chains[[4]] <- map(pRvsPr, function(x) MCMCchains(x$mcmc_samples, params = "beta"))


beta_prob <- map(beta_chains, function(x) map_dbl(x, function(y) length(which(y>0))/length(y)))

df1 <- data.frame(beta = unlist(beta_prob),
                  type = c(rep("R", length(beta_chains[[1]])),
                           rep("N", length(beta_chains[[2]])),
                           rep("highN", length(beta_chains[[3]])),
                           rep("PR", length(beta_chains[[4]]))))

theme_set(theme_bw())
g1 <- ggplot(df1, aes(x = type, y = beta)) + 
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="orange") +
  scale_x_discrete("type", labels = c("N (95%)", 
                                      bquote(bar("N")), 
                                      bquote("P("*bar("r")*">0)"),
                                      bquote(bar("r")))) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,1)) +
  labs(y = bquote("P("*beta*">0)")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


R_sq <- list()
R_sq[[1]] <- map_dbl(RvsPr, function(x) x$mcmc_sum[1,1])
R_sq[[2]] <- map_dbl(NvsPr, function(x) x$mcmc_sum[1,1])
R_sq[[3]] <- map_dbl(hNvsPr, function(x) x$mcmc_sum[1,1])
R_sq[[4]] <- map_dbl(pRvsPr, function(x) x$mcmc_sum[1,1])

df2 <- data.frame(R_sq = unlist(R_sq),
                  type = c(rep("R", length(R_sq[[1]])),
                           rep("N", length(R_sq[[2]])),
                           rep("highN", length(R_sq[[3]])),
                           rep("PR", length(R_sq[[4]]))))

g2 <- ggplot(df2, aes(x = type, y = R_sq)) + 
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5, 
               fill="orange") +
  scale_x_discrete("type", labels = c("N (97,5%)", 
                                      bquote(bar("N")), 
                                      bquote("P("*bar("r")*">0)"),
                                      bquote(bar("r")))) +
  scale_y_continuous(breaks=seq(0, 0.3, 0.05), limits = c(0,0.32)) +
  labs(y = bquote("R"^2)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

g1 + g2 + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("paper_results/fig3.tiff", width = 8, height = 6, units = "in")


# Discrimination ability of SDMs ------------------------------------------

# Pocc of maps sites
ss <- map_dbl(results_maps_hmsc, length)
species <- map2(spcode, ss, function(x,y) rep(x, times = y)) %>% unlist()

df3 <- data.frame(poc = unlist(results_maps_hmsc),
                  species = as.factor(species))

g3 <- ggplot(df3, aes(x = species, y = poc)) + 
  geom_boxplot(fill = "grey", alpha = 0.5, col = "black") +
  labs(y = bquote("P"[occ])) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 8, angle = 90))

# AUC and Pocc comparison
medPr <- map_dbl(results_maps_hmsc, median)

g4 <- ggplot(mapping = aes(x = results_auc, y = medPr)) +
  geom_point(size = 2.5, color = "orange") +
  labs(x = "AUC", y = bquote("Median" ~ "P"[occ])) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  scale_x_continuous(breaks=seq(0.65, 1, 0.05), limits = c(0.65, 1)) +
  scale_y_continuous(breaks=seq(0.50, 1, 0.05), limits = c(0.55, 1))

# Discrimination between rmax<0 and rmax>0
pa <- map(avgR, function(x) ifelse(x>=1,"r>0","r<0"))
n_pa <- map_dbl(pa, function(x) c(length(which(x=="r<0"))))
index <- which(n_pa>=5)

poc2 <- results_maps_hmsc[index]
pa2 <- pa[index]
ss2 <- ss[index]
spcode2 <- spcode[index]
species2 <- map2(spcode2, ss2, function(x,y) rep(x, times = y)) %>% unlist()
n_pa2 <- map(pa2, function(x) c(length(which(x=="r<0")), length(which(x=="r>0"))))

df4 <- data.frame(poc = unlist(poc2),
                  pa = as.factor(unlist(pa2)),
                  species = as.factor(species2))

df_text <- data.frame(n = unlist(n_pa2),
                      species = rep(spcode2, each = 2),
                      r = rep(c("r<0", "r>0"), times = length(n_pa2)))

g5 <- ggplot(df4, aes(x = species, y = poc)) + 
  geom_boxplot(aes(fill = pa), alpha = 0.5, col = "black") +
  labs(y = bquote("P"[occ])) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  scale_fill_manual(name = "Growth Rate", values = c("r<0" = "grey", "r>0" = "orange")) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,1.03)) +
  geom_text(data = df_text, aes(x = species, y= 1.03, label = n), 
            position=position_dodge2(width=0.7),
            size = 2)


# Load packages and data --------------------------------------------------

library(tidyverse)
library(MCMCvis)
library(doSNOW)
library(foreach)
library(rstan)
library(MCMCvis)
library(ggthemes)
library(boot)
library(grid)
library(sf)

# Install with devtools::install_github("thomasp85/patchwork")
library(patchwork)

results_hmsc_maps <- readRDS("results_hmsc_maps.rds")
results_brt_maps <- readRDS("results_brt_maps.rds")
results_cjspop <- readRDS("results_cjspop.rds")
results_dem <- readRDS("results_dem.rds")
results_pR <- readRDS("results_pR.rds")
auc_hmsc <- readRDS("results_auc_hmsc.rds")
auc_brt <- readRDS("results_auc_brt.rds")
data_sdm <- readRDS("data_sdm_pca.rds")
occ <- readRDS("data_occ.rds")

spcode <- results_cjspop$spcode

N <- results_dem$N
pR <- list()
for (i in 1:length(results_pR)) {
  
  pR[[i]] <- results_pR[[i]]
  if (any(pR[[i]]==1)) pR[[i]][which(pR[[i]]==1)] <- 0.9999
  
}

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

x_std <- map(results_hmsc_maps, function(x) (x - mean(x))/sd(x))

## R vs Pr
RvsPr <- foreach (i=1:length(popR)) %do% {
  
  stan(file = 'lm.stan', 
       data = list(X = x_std[[i]], 
                   y = log(popR[[i]]), 
                   N = length(popR[[i]])),
       iter = 2000)
  
}

## N vs Pr
NvsPr <- foreach (i=1:length(popR)) %do% {
  
  stan(file = 'lm.stan', 
       data = list(X = x_std[[i]], 
                   y = log(N[[i]]$popN), 
                   N = length(popR[[i]])),
       iter = 2000)
  
}

## High N vs Pr
hNvsPr <- foreach (i=1:length(popR)) %do% {
  
  stan(file = 'lm.stan', 
       data = list(X = x_std[[i]], 
                   y = log(N[[i]]$highN), 
                   N = length(popR[[i]])),
       iter = 2000)
  
}


## PR vs Pr
pRvsPr <- foreach (i=1:(length(popR))) %do% {
  
  stan(file = 'lm_beta.stan', 
       data = list(X = x_std[[i]], 
                   y = pR[[i]], 
                   N = length(popR[[i]])),
       iter = 2000)
  
}

## Check convergence
map_dbl(RvsPr, function(x) max(MCMCsummary(x)[,"Rhat"]))
map_dbl(NvsPr, function(x) max(MCMCsummary(x)[,"Rhat"]))
map_dbl(hNvsPr, function(x) max(MCMCsummary(x)[,"Rhat"]))
map_dbl(pRvsPr, function(x) max(MCMCsummary(x)[,"Rhat"]))

map_dbl(RvsPr, function(x) min(MCMCsummary(x)[,"n.eff"]))
map_dbl(NvsPr, function(x) min(MCMCsummary(x)[,"n.eff"]))
map_dbl(hNvsPr, function(x) min(MCMCsummary(x)[,"n.eff"]))
map_dbl(pRvsPr, function(x) min(MCMCsummary(x)[,"n.eff"]))


## Plot
beta_chains <-  list()
beta_chains[[1]] <- map(RvsPr, function(x) MCMCchains(x, params = "beta"))
beta_chains[[2]] <- map(NvsPr, function(x) MCMCchains(x, params = "beta"))
beta_chains[[3]] <- map(hNvsPr, function(x) MCMCchains(x, params = "beta"))
beta_chains[[4]] <- map(pRvsPr, function(x) MCMCchains(x, params = "beta"))

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
               dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c("N (95%)", 
                                      bquote(bar("N")), 
                                      bquote("P("*bar("r")*">0)"),
                                      bquote(bar("r")))) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,1)) +
  labs(y = bquote("P("*beta*">0)"), title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


R_sq <- list()
R_sq[[1]] <- map_dbl(RvsPr, function(x) MCMCsummary(x, params = "Rsq")[1,1])
R_sq[[2]] <- map_dbl(NvsPr, function(x) MCMCsummary(x, params = "Rsq")[1,1])
R_sq[[3]] <- map_dbl(hNvsPr, function(x) MCMCsummary(x, params = "Rsq")[1,1])
R_sq[[4]] <- map_dbl(pRvsPr, function(x)MCMCsummary(x, params = "Rsq")[1,1])

df2 <- data.frame(R_sq = unlist(R_sq),
                  type = c(rep("R", length(R_sq[[1]])),
                           rep("N", length(R_sq[[2]])),
                           rep("highN", length(R_sq[[3]])),
                           rep("PR", length(R_sq[[4]]))))

g2 <- ggplot(df2, aes(x = type, y = R_sq)) + 
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c("N (95%)", 
                                      bquote(bar("N")), 
                                      bquote("P("*bar("r")*">0)"),
                                      bquote(bar("r")))) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,0.7)) +
  labs(y = bquote("R"^2), title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


# Repeat analysis with BRT results ----------------------------------------

x_std_brt <- map(results_brt_maps, function(x) (x - mean(x))/sd(x))

## R vs Pr
RvsPr_brt <- foreach (i=1:length(popR)) %do% {
  
  stan(file = 'lm.stan', 
       data = list(X = x_std_brt[[i]], 
                   y = log(popR[[i]]), 
                   N = length(popR[[i]])),
       iter = 2000)
  
}

## N vs Pr
NvsPr_brt <- foreach (i=1:length(popR)) %do% {
  
  stan(file = 'lm.stan', 
       data = list(X = x_std_brt[[i]], 
                   y = log(N[[i]]$popN), 
                   N = length(popR[[i]])),
       iter = 2000)
  
}

## High N vs Pr
hNvsPr_brt <- foreach (i=1:length(popR)) %do% {
  
  stan(file = 'lm.stan', 
       data = list(X = x_std_brt[[i]], 
                   y = log(N[[i]]$highN), 
                   N = length(popR[[i]])),
       iter = 2000)
  
}


## PR vs Pr
pRvsPr_brt <- foreach (i=1:(length(popR))) %do% {
  
  stan(file = 'lm_beta.stan', 
       data = list(X = x_std_brt[[i]], 
                   y = pR[[i]], 
                   N = length(popR[[i]])),
       iter = 2000)
  
}

## Check convergence
map_dbl(RvsPr_brt, function(x) max(MCMCsummary(x)[,"Rhat"]))
map_dbl(NvsPr_brt, function(x) max(MCMCsummary(x)[,"Rhat"]))
map_dbl(hNvsPr_brt, function(x) max(MCMCsummary(x)[,"Rhat"]))
map_dbl(pRvsPr_brt, function(x) max(MCMCsummary(x)[,"Rhat"]))

map_dbl(RvsPr_brt, function(x) min(MCMCsummary(x)[,"n.eff"]))
map_dbl(NvsPr_brt, function(x) min(MCMCsummary(x)[,"n.eff"]))
map_dbl(hNvsPr_brt, function(x) min(MCMCsummary(x)[,"n.eff"]))
map_dbl(pRvsPr_brt, function(x) min(MCMCsummary(x)[,"n.eff"]))


## Plot

beta_chains_brt <-  list()
beta_chains_brt[[1]] <- map(RvsPr_brt, function(x) MCMCchains(x, params = "beta"))
beta_chains_brt[[2]] <- map(NvsPr_brt, function(x) MCMCchains(x, params = "beta"))
beta_chains_brt[[3]] <- map(hNvsPr_brt, function(x) MCMCchains(x, params = "beta"))
beta_chains_brt[[4]] <- map(pRvsPr_brt, function(x) MCMCchains(x, params = "beta"))

beta_prob_brt <- map(beta_chains_brt, function(x) map_dbl(x, function(y) length(which(y>0))/length(y)))

df3 <- data.frame(beta = unlist(beta_prob_brt),
                  type = c(rep("R", length(beta_chains_brt[[1]])),
                           rep("N", length(beta_chains_brt[[2]])),
                           rep("highN", length(beta_chains_brt[[3]])),
                           rep("PR", length(beta_chains_brt[[4]]))))

theme_set(theme_bw())
g3 <- ggplot(df3, aes(x = type, y = beta)) + 
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .4,
               fill="orange") +
  scale_x_discrete("type", labels = c("N (95%)", 
                                      bquote(bar("N")), 
                                      bquote("P("*bar("r")*">0)"),
                                      bquote(bar("r")))) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,1)) +
  labs(y = bquote("P("*beta*">0)"), title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


R_sq_brt <- list()
R_sq_brt[[1]] <- map_dbl(RvsPr_brt, function(x) MCMCsummary(x, params = "Rsq")[1,1])
R_sq_brt[[2]] <- map_dbl(NvsPr_brt, function(x) MCMCsummary(x, params = "Rsq")[1,1])
R_sq_brt[[3]] <- map_dbl(hNvsPr_brt, function(x) MCMCsummary(x, params = "Rsq")[1,1])
R_sq_brt[[4]] <- map_dbl(pRvsPr_brt, function(x)MCMCsummary(x, params = "Rsq")[1,1])

df4 <- data.frame(R_sq = unlist(R_sq_brt),
                  type = c(rep("R", length(R_sq_brt[[1]])),
                           rep("N", length(R_sq_brt[[2]])),
                           rep("highN", length(R_sq_brt[[3]])),
                           rep("PR", length(R_sq_brt[[4]]))))

g4 <- ggplot(df4, aes(x = type, y = R_sq)) + 
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c("N (95%)", 
                                      bquote(bar("N")), 
                                      bquote("P("*bar("r")*">0)"),
                                      bquote(bar("r")))) +
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,0.7)) +
  labs(y = bquote("R"^2), title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


g1 + g3 + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("paper_results2/sdm_fig1.jpeg", width = 8, height = 6, units = "in")

g2 + g4 + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("paper_results2/sdm_fig2.jpeg", width = 8, height = 6, units = "in")

# Comparison of Pocc with r<0 and r>0 -------------------------------------

# Select species with negative r populations
neg_r <- map_dbl(popR, function(x) sum(log(x)<0))
pos_r <- map_dbl(popR, function(x) sum(log(x)>0))

ss <- apply(cbind(neg_r, pos_r), 1, function(x) paste(x[1], x[2], sep=" - ")) 

popR2 <- popR[neg_r>=4]
pocc_hmsc <- results_hmsc_maps[neg_r>=4]
pocc_brt <- results_brt_maps[neg_r>=4]

spcode2 <- spcode[neg_r>=4]
spcode3 <- map2(spcode2, map_dbl(popR2, length), function(x,y) rep(x, y)) %>%
  unlist()

r_group <- map(popR2, function(x) as.factor(as.numeric(log(x)>0)))

p_brt <- c()
p_hmsc <- c()
for (i in 1:length(pocc_brt)) {
  
  p_brt[i] <- t.test(pocc_brt[[i]]~r_group[[i]], alternative = "less")$p.value
  p_hmsc[i] <- t.test(pocc_hmsc[[i]]~r_group[[i]], alternative = "less")$p.value
  
}

box_data1 <- data.frame(y = unlist(pocc_hmsc),
                        x = spcode3,
                        group = unlist(r_group))

box_data2 <- data.frame(y = unlist(pocc_brt),
                        x = spcode3,
                        group = unlist(r_group))

g5 <- ggplot() +
  geom_boxplot(mapping = aes(y = y, x = x, fill = group), data = box_data1, alpha = 0.7) +
  geom_text(mapping = aes(x = spcode2, y = 1.03, label = round(p_hmsc, 3)), size = 3) +
  geom_text(mapping = aes(x = spcode2, y = -0.05, label = ss[neg_r>=4]), size = 3) +
  scale_fill_manual(labels = c(bquote(bar("r")*"<0"), bquote(bar("r")*">0")), 
                    values = c("grey", "orange")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  labs(y = bquote("P"[occ]), title = "HMSC")

g6 <- ggplot() +
  geom_boxplot(mapping = aes(y = y, x = x, fill = group), data = box_data2, alpha = 0.7) +
  geom_text(mapping = aes(x = spcode2, y = 1.03, label = round(p_brt, 3)), size = 3) +
  geom_text(mapping = aes(x = spcode2, y = -0.05, label = ss[neg_r>=4]), size = 3) +
  scale_fill_manual(labels = c(bquote(bar("r")*"<0"), bquote(bar("r")*">0")), 
                    values = c("grey", "orange")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  labs(y = bquote("P"[occ]), title = "BRT")

g5/g6
ggsave("paper_results2/sdm_fig3.jpeg", width = 8, height = 8, units = "in")


# Other Supplementary Figures ---------------------------------------------

beta_r_hmsc <- map_dbl(RvsPr, function(x) x$mcmc_sum["beta",1])
beta_r_brt <-  map_dbl(RvsPr_brt, function(x) x$mcmc_sum["beta",1])

sp_eng2 <- sp_eng
sp_eng2[which(sp_eng == "Brown Creeper")] <- "American Treecreeper"

traits <- read.delim("BirdFuncDat.txt", na.strings = " ") %>%
  as.data.frame() %>%
  filter(English %in% sp_eng2)


diets <- traits$Diet.5Cat[order(traits$English)]
bdsize <- traits$BodyMass.Value[order(traits$English)]

beta_r_hmsc2 <- beta_r_hmsc[order(sp_eng2)]
beta_r_brt2 <- beta_r_brt[order(sp_eng2)]

r_sign_hmsc2 <- c()
r_sign_brt2 <- c()
for (i in 1:length(beta_r_hmsc)) {
  
  r_sign_hmsc2[i] <- ifelse(beta_r_hmsc2[i]<0, 1, 0)
  r_sign_brt2[i] <- ifelse(beta_r_brt2[i]<0, 1, 0)
  
}

sup_df2_hmsc <- data.frame(beta_r = beta_r_hmsc2,
                           diets = diets,
                           bdsize = bdsize,
                           sign_r = as.factor(r_sign_hmsc2))

sup_df2_brt <- data.frame(beta_r = beta_r_brt2,
                           diets = diets,
                           bdsize = bdsize,
                           sign_r = as.factor(r_sign_brt2))

## Diet

sg3 <- ggplot(data = sup_df2_hmsc, aes(x = diets, y = beta_r)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
                                      #bquote(beta*"<0"))) +
  labs(y = bquote(beta), title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

sg4 <- ggplot(data = sup_df2_brt, aes(x = diets, y = beta_r)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
  #bquote(beta*"<0"))) +
  labs(y = bquote(beta), title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

## Body size

sg5 <- ggplot(data = sup_df2_hmsc, aes(x = sign_r, y = bdsize)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "Body Size", title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

sg6 <- ggplot(data = sup_df2_brt, aes(x = sign_r, y = bdsize)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "Body Size", title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

## Nesting Behaviour

nests <- c("ground", "tree", "shrub", "cavity", "tree", "tree", "cavity", "cavity",
           "tree", "ground", "cavity", "shrub", "cavity", "tree", "shrub", "ground",
           "shrub")

sup_df3_hmsc <- data.frame(beta_r = beta_r_hmsc,
                           nests = nests)

sup_df3_brt <- data.frame(beta_r = beta_r_brt,
                           nests = nests)

sg7 <- ggplot(data = sup_df3_hmsc, aes(x = nests, y = beta_r)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
  #bquote(beta*"<0"))) +
  labs(y = bquote(beta), title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

sg8 <- ggplot(data = sup_df3_brt, aes(x = nests, y = beta_r)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
  #bquote(beta*"<0"))) +
  labs(y = bquote(beta), title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

(sg3 + sg4) / (sg7 + sg8) / (sg5 + sg6) + 
  plot_annotation(tag_levels = 'a')
ggsave("paper_results2/sup_fig1.jpeg", width = 8, height = 10, units = "in")

## sample size and AUC

nch <- map_dbl(chdata, function(x) nrow(x$ch_year))
npa <- map_dbl(data_sdm, function(x) nrow(x)/2)

sup_df4_hmsc <- data.frame(nch = nch,
                           npa = npa,
                           auc = auc_hmsc,
                           sign_r = as.factor(r_sign_hmsc))

sup_df4_brt <- data.frame(nch = nch,
                          npa = npa,
                          auc = auc_brt,
                          sign_r = as.factor(r_sign_brt))

sg9 <- ggplot(data = sup_df4_hmsc, aes(x = sign_r, y = nch)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "Number of Individuals", title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

sg10 <- ggplot(data = sup_df4_brt, aes(x = sign_r, y = nch)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "Number of Individuals", title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

sg11 <- ggplot(data = sup_df4_hmsc, aes(x = sign_r, y = npa)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "Number of Presences", title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

sg12 <- ggplot(data = sup_df4_brt, aes(x = sign_r, y = npa)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "Number of Presences", title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

sg13 <- ggplot(data = sup_df4_hmsc, aes(x = sign_r, y = auc)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "AUC", title = "HMSC") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

sg14 <- ggplot(data = sup_df4_brt, aes(x = sign_r, y = auc)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .4, 
               fill="orange") +
  scale_x_discrete("type", labels = c(bquote(beta*">0"), 
                                      bquote(beta*"<0"))) +
  labs(y = "AUC", title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 10))

(sg9 + sg10) / (sg11 + sg12) / (sg13 + sg14) + 
  plot_annotation(tag_levels = 'a')
ggsave("paper_results2/sup_fig2.jpeg", width = 8, height = 10, units = "in")


# Tables ------------------------------------------------------------------

idx <- order(spcode)
np <- map_dbl(data_sdm, function(x) nrow(x)/2)
ab_range <- c("US", "West", "West", "US", "US", "West", "West", "East",
              "East", "US", "East", "West", "US", "US", "US", "US", "US")

table1 <- data.frame(np = np[idx],
                     range = ab_range[idx],
                     auc_hmsc = round(auc_hmsc[idx],2),
                     auc_brt = round(auc_brt[idx],2))

write.csv(table1, file = "paper_results2/table1a.csv")



table2a <- data.frame(spcode = spcode[idx],
                      prob_hmsc = round(beta_prob[[1]][idx], 3),
                      prob_brt = round(beta_prob_brt[[1]][idx], 3),
                      Rsq_hmsc = round(R_sq[[1]][idx], 2),
                      Rsq_brt = round(R_sq_brt[[1]][idx], 2))

write.csv(table2a, file = "paper_results2/table2a.csv")

table2b <- data.frame(spcode = spcode[idx],
                      prob_hmsc = round(beta_prob[[2]][idx], 3),
                      prob_brt = round(beta_prob_brt[[2]][idx], 3),
                      Rsq_hmsc = round(R_sq[[2]][idx], 2),
                      Rsq_brt = round(R_sq_brt[[2]][idx], 2))

write.csv(table2b, file = "paper_results2/table2b.csv")




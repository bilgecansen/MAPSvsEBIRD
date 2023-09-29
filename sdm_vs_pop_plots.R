
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
library(pROC)
library(raster)
library(factoextra)
library(gbm)

# Install with devtools::install_github("thomasp85/patchwork")
library(patchwork)

results_brt <- readRDS("results_brt.rds")
results_hmsc_maps <- readRDS("results_hmsc_maps.rds")
results_brt_maps <- readRDS("results_brt_maps.rds")
results_cjspop <- readRDS("results_cjspop.rds")
results_dem <- readRDS("results_dem.rds")
results_pR <- readRDS("results_pR.rds")
auc_hmsc_full <- readRDS("results_auc_hmsc.rds")
auc_brt_full <- readRDS("results_auc_brt.rds")
data_sdm <- readRDS("data_sdm_pca.rds")
occ <- readRDS("data_occ.rds")
chdata <- results_cjspop$chdata

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
  labs(y = bquote("P("*beta*">0)"), title = "GLM") +
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
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,0.45)) +
  labs(y = bquote("R"^2), title = "GLM") +
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
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0,0.45)) +
  labs(y = bquote("R"^2), title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


g1 + g3 + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("fig1.jpeg", width = 16.8, height = 12, units = "cm")

g2 + g4 + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("fig2.jpeg", width = 16.8, height = 12, units = "cm")


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

auc_brt <- map2_dbl(r_group, pocc_brt, function(x,y) auc(x, y, direction = "<"))
auc_hmsc <- map2_dbl(r_group, pocc_hmsc, function(x,y) auc(x, y, direction = "<"))

auc_random_brt <- map2(r_group, pocc_brt, function(x,y) replicate(1000, auc(sample(x), y, direction = "<")))
p_value_brt <- map2_dbl(auc_random_brt, auc_brt, function(x,y) length(which(x > y))/1000)

auc_random_hmsc <- map2(r_group, pocc_hmsc, function(x,y) replicate(1000, auc(sample(x), y, direction = "<")))
p_value_hmsc <- map2_dbl(auc_random_hmsc, auc_hmsc, function(x,y) length(which(x > y))/1000)

auc_brt_text <- paste(round(auc_brt, 3), map_chr(p_value_brt, function(x) ifelse(x <0.05, "*", "")), sep = "")
auc_hmsc_text <- paste(round(auc_hmsc, 3), map_chr(p_value_hmsc, function(x) ifelse(x <0.05, "*", "")), sep = "")

box_data1 <- data.frame(y = unlist(pocc_hmsc),
                        x = spcode3,
                        group = unlist(r_group))

box_data2 <- data.frame(y = unlist(pocc_brt),
                        x = spcode3,
                        group = unlist(r_group))

g5 <- ggplot() +
  geom_boxplot(mapping = aes(y = y, x = x, fill = group), data = box_data1, alpha = 0.7) +
  geom_text(mapping = aes(x = spcode2, y = 1.03, label = auc_hmsc_text), size = 3) +
  geom_text(mapping = aes(x = spcode2, y = -0.05, label = ss[neg_r>=4]), size = 3) +
  scale_fill_manual(labels = c(bquote(bar("r")*"<0"), bquote(bar("r")*">0")), 
                    values = c("grey", "orange")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Suitability", title = "GLM")

g6 <- ggplot() +
  geom_boxplot(mapping = aes(y = y, x = x, fill = group), data = box_data2, alpha = 0.7) +
  geom_text(mapping = aes(x = spcode2, y = 1.03, label = auc_brt_text), size = 3) +
  geom_text(mapping = aes(x = spcode2, y = -0.05, label = ss[neg_r>=4]), size = 3) +
  scale_fill_manual(labels = c(bquote(bar("r")*"<0"), bquote(bar("r")*">0")), 
                    values = c("grey", "orange")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Suitability", title = "BRT")

g5/g6 +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("fig3.jpeg", width = 16.8, height = 16.8, units = "cm")


# Maps --------------------------------------------------------------------

# WE don't have permission to share band and raster data via a repository
# please contact bilgecan.sen@gmail.com if you need access to this data.

data_band <- readRDS("band_data.rds")
rasters <- paste("weather/average", list.files("weather/average"), sep = "/")
r_avg <- stack(rasters)
r_avg <- rotate(r_avg)

z <- values(r_avg)
coord <- coordinates(r_avg)

plot_map <- function(i, xlimits, ylimits, title, lpos, f = NULL) {
  eig <- get_eigenvalue(chdata[[i]]$res_pca)
  dim_num <- which(eig$cumulative.variance.percent>80)[1]
  z_pca <- predict(chdata[[i]]$res_pca, z)[,1:dim_num]
  colnames(z_pca) <- paste("PC", 1:dim_num, sep = "")
  idx_na <- which(is.na(z_pca[,1]))
  
  # Pop coordinates
  pop <- dplyr::select(data_band, pop, lat, long) %>%
    distinct(pop, lat, long) %>%
    filter(lat < 99 | long < 0) %>%
    group_by(pop) %>%
    summarise(lat = mean(lat),
              long = mean(long))
  
  pop_coord <- filter(pop, pop %in% as.numeric(row.names(chdata[[i]]$Nobs_ad)))
  idx_long <- pop_coord$long < 0
  pop_coord <- pop_coord[idx_long,]
  pop_coord$r <- popR[[i]][idx_long]
 
  pop_coord$r <- ifelse(pop_coord$r < 1, "sink", "source")
  
  # BRT predictions
  z_pred <- predict.gbm(results_brt[[i]], as.data.frame(z_pca)[-idx_na,], type = "response")
  
  theme_set(theme_bw())
  ggplot() +
    geom_raster(aes(x = coord[-idx_na,1], y = coord[-idx_na,2], fill = z_pred)) +
    coord_quickmap() +
    scale_fill_gradient(limits = c(0,1),
                        low = "grey", 
                        high = "green4") +
    geom_point(data = pop_coord, mapping = aes(x = long, y = lat, col = r), size = 2, alpha = 0.7) +
    scale_color_manual(name = NULL, 
                       values = c("sink" = "darkred", 
                                  "source" = "blue4")) +
    scale_x_continuous(limits = xlimits) +
    scale_y_continuous(limits = ylimits) +
    labs(title = title) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.key = element_rect(fill = "gray"),
          legend.position = lpos)
}

g7 <- plot_map(16, xlimits = c(-125, -65), ylimits = c(25, 50), 
               title = "Spotted Towhee (SPTO)", lpos = "none")

g8 <- plot_map(4, xlimits = c(-125, -65), ylimits = c(25, 50), title = "Hairy Woodpecker (HAWO)", lpos = "bottom")

g7/g8 +
  plot_annotation(tag_levels = 'a')
ggsave("fig4.jpeg", width = 16.8, height = 23, units = "cm")

plot_map(1, xlimits = c(-125, -65), ylimits = c(25, 50), title = "RCSP", lpos = "bottom")
ggsave("fig_sup_RCSP.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(2, xlimits = c(-125, -65), ylimits = c(25, 50), title = "HUVI", lpos = "bottom")
ggsave("fig_sup_HUVI.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(3, xlimits = c(-125, -65), ylimits = c(25, 50), title = "CALT", lpos = "bottom")
ggsave("fig_sup_CALT.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(5, xlimits = c(-125, -65), ylimits = c(25, 50), title = "BRCR", lpos = "bottom")
ggsave("fig_sup_BRCR.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(6, xlimits = c(-125, -65), ylimits = c(25, 50), title = "WEWP", lpos = "bottom")
ggsave("fig_sup_WEWP.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(7, xlimits = c(-125, -65), ylimits = c(25, 50), title = "MOCH", lpos = "bottom")
ggsave("fig_sup_MOCH.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(8, xlimits = c(-125, -65), ylimits = c(25, 50), title = "CACH", lpos = "bottom")
ggsave("fig_sup_CACH.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(14, xlimits = c(-125, -65), ylimits = c(25, 50), title = "REVI", lpos = "bottom")
ggsave("fig_sup_REVI.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)

plot_map(15, xlimits = c(-125, -65), ylimits = c(25, 50), title = "INBU", lpos = "bottom")
ggsave("fig_sup_INBU.pdf", width = 16.8, height = 16.8, units = "cm", dpi = 600)


# Other Supplementary Figures ---------------------------------------------

spname <- c("Aimophila ruficeps", "Vireo huttoni","Pipilo crissalis",
            "Picoides villosus", "Certhia americana", "Contopus sordidulus",
            "Parus gambeli", "Parus carolinensis", "Empidonax virescens",
            "Melospiza lincolnii", "Baeolophus bicolor", "Chamaea fasciata",
            "Picoides pubescens", "Vireo olivaceus", "Passerina cyanea",
            "Pipilo maculatus", "Icteria virens")

traits <- read.delim("BirdFuncDat.txt", na.strings = " ") %>%
  as.data.frame() %>%
  filter(Scientific %in% spname2)

diets <- traits$Diet.5Cat[order(traits$Scientific)]
bdsize <- traits$BodyMass.Value[order(traits$Scientific)]

sup_df2_hmsc <- data.frame(auc = auc_hmsc,
                           diets = diets,
                           bdsize = bdsize)

sup_df2_brt <- data.frame(auc = auc_brt,
                          diets = diets,
                          bdsize = bdsize)

nests2 <- nests[neg_r>=4]

## Diet
sg3 <- ggplot(data = sup_df2_hmsc, aes(x = diets, y = auc)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
                                      #bquote(beta*"<0"))) +
  labs(y = "AUC (Source vs Sink)", title = "GLM") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

sg4 <- ggplot(data = sup_df2_brt, aes(x = diets, y = auc)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
  #bquote(beta*"<0"))) +
  labs(y = "AUC (Source vs Sink)", title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))

## Body size
stan(file = 'lm_beta.stan', 
     data = list(X = ((sup_df2_hmsc$bdsize - mean(sup_df2_hmsc$bdsize))/sd(sup_df2_hmsc$bdsize)), 
                 y = sup_df2_hmsc$auc, 
                 N = length(sup_df2_hmsc$auc)),
     iter = 2000)
sg5 <- ggplot(data = sup_df2_hmsc, aes(x = bdsize, y = auc)) +
  geom_point() +
  geom_smooth(se = T, method = "lm") +
  labs(x = "Body Size (g)", y = "AUC (Source vs Sink)", 
       title = "GLM (R-squared = 0.21)") +
  theme(axis.text.x = element_text(size = 8))


stan(file = 'lm_beta.stan', 
     data = list(X = ((sup_df2_brt$bdsize - mean(sup_df2_brt$bdsize))/sd(sup_df2_brt$bdsize)), 
                 y = sup_df2_brt$auc, 
                 N = length(sup_df2_brt$auc)),
     iter = 2000)
sg6 <- ggplot(data = sup_df2_brt, aes(x = bdsize, y = auc)) +
  geom_point() +
  geom_smooth(se = T, method = "lm") +
  labs(x = "Body Size (g)", y = "AUC (Source vs Sink)",  
       title = "BRT (R-squared = 0.21)") +
  theme(axis.text.x = element_text(size = 8))

## Nesting Behaviour
nests <- c("ground", "tree", "shrub", "cavity", "tree", "tree", "cavity", "cavity",
           "tree", "ground", "cavity", "shrub", "cavity", "tree", "shrub", "ground",
           "shrub")

sup_df3_hmsc <- data.frame(auc = auc_hmsc,
                           nests = nests2)

sup_df3_brt <- data.frame(auc = auc_brt,
                          nests = nests2)

sg7 <- ggplot(data = sup_df3_hmsc, aes(x = nests, y = auc)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
  #bquote(beta*"<0"))) +
  labs(y = "AUC (Source vs Sink)", title = "GLM") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10))

sg8 <- ggplot(data = sup_df3_brt, aes(x = nests, y = auc)) +
  geom_boxplot(fill = "grey", alpha = 0.5) +
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               #dotsize = .5, 
               fill="orange") +
  #scale_x_discrete("type", labels = c(bquote(beta*">0"),
  #bquote(beta*"<0"))) +
  labs(y = "AUC (Source vs Sink)", title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10))

(sg3 + sg4) / (sg7 + sg8) / (sg5 + sg6) + 
  plot_annotation(tag_levels = 'a')
ggsave("sup_fig1.jpeg", width = 8, height = 10, units = "in")

## sample size and AUC
nch <- map_dbl(chdata, function(x) nrow(x$ch_year))
npa <- map_dbl(data_sdm, function(x) nrow(x)/2)

sup_df4_hmsc <- data.frame(nch = nch[neg_r>=4],
                           npa = npa[neg_r>=4],
                           auc_full = auc_hmsc_full[neg_r>=4],
                           auc = auc_hmsc)

sup_df4_brt <- data.frame(nch = nch[neg_r>=4],
                          npa = npa[neg_r>=4],
                          auc_full = auc_brt_full[neg_r>=4],
                          auc = auc_brt)


stan(file = 'lm_beta.stan', 
     data = list(X = ((nch - mean(nch))/sd(nch))[neg_r>=4], 
                 y = sup_df4_hmsc$auc, 
                 N = length(sup_df4_hmsc$auc)),
     iter = 2000)
sg9 <- ggplot(data = sup_df4_hmsc, aes(y = auc, x = nch)) +
  geom_point() +
  geom_smooth(se = T, method = "lm") +
  labs(x = "Number of captured individuals", y = "AUC (Source vs Sink)",
       title = "GLM (R-squared = 0.15)") +
  theme(axis.text.x = element_text(size = 8))

stan(file = 'lm_beta.stan', 
     data = list(X = ((nch - mean(nch))/sd(nch))[neg_r>=4], 
                 y = sup_df4_brt$auc, 
                 N = length(sup_df4_brt$auc)),
     iter = 2000)
sg10 <- ggplot(data = sup_df4_brt, aes(y = auc, x = nch)) +
  geom_point() +
  geom_smooth(se = T, method = "lm") +
  labs(x = "Number of captured individuals", y = "AUC (Source vs Sink)",
       title = "BRT (R-squared = 0.39)") +
  theme(axis.text.x = element_text(size = 8))

stan(file = 'lm_beta.stan', 
     data = list(X = ((npa - mean(npa))/sd(npa))[neg_r>=4], 
                 y = sup_df4_hmsc$auc, 
                 N = length(sup_df4_hmsc$auc)),
     iter = 2000)
sg11 <- ggplot(data = sup_df4_hmsc, aes(y = auc, x = npa)) +
  geom_point() +
  geom_smooth(se = T, method = "lm") +
  labs(x = "Number of presences", y = "AUC (Source vs Sink)",
       title = "GLM (R-squared = 0.15)") +
  theme(axis.text.x = element_text(size = 8))

stan(file = 'lm_beta.stan', 
     data = list(X = ((npa - mean(npa))/sd(npa))[neg_r>=4], 
                 y = sup_df4_brt$auc, 
                 N = length(sup_df4_brt$auc)),
     iter = 2000)
sg12 <- ggplot(data = sup_df4_brt, aes(y = auc, x = npa)) +
  geom_point() +
  geom_smooth(se = T, method = "lm") +
  labs(x = "Number of presences", y = "AUC (Source vs Sink)",
       title = "BRT (R-squared = 0.18)") +
  theme(axis.text.x = element_text(size = 8))

x <- ((auc_hmsc_full - mean(auc_hmsc_full))/sd(auc_hmsc_full))[neg_r>=4]
stan(file = 'lm_beta2.stan', 
     data = list(X = cbind(x, x^2), 
                 y = sup_df4_hmsc$auc, 
                 N = length(sup_df4_hmsc$auc)),
     iter = 2000)
sg13 <- ggplot(data = sup_df4_hmsc, aes(y = auc, x = auc_full)) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2)) +
  labs(x = "AUC (Presence vs Pseudoabsence)", y = "AUC (Source vs Sink)",
       title = "GLM (R-squared = 0.33)") +
  theme(axis.text.x = element_text(size = 8))

x <- ((auc_brt_full - mean(auc_brt_full))/sd(auc_brt_full))[neg_r>=4]
stan(file = 'lm_beta2.stan', 
     data = list(X = cbind(x, x^2), 
                 y = sup_df4_brt$auc, 
                 N = length(sup_df4_brt$auc)),
     iter = 2000)
sg14 <- ggplot(data = sup_df4_brt, aes(y = auc, x = auc_full)) +
  geom_point() +
  geom_smooth(se = T, method = "lm", formula = y ~ poly(x, 2)) +
  labs(x = "AUC (Presence vs Pseudoabsence)", y = "AUC (Source vs Sink)",
       title = "BRT (R-squared = 0.48)") +
  theme(axis.text.x = element_text(size = 8))

(sg9 + sg10) / (sg11 + sg12) / (sg13 + sg14) + 
  plot_annotation(tag_levels = 'a')
ggsave("sup_fig2.jpeg", width = 8, height = 10, units = "in")


# Tables ------------------------------------------------------------------

idx <- order(spcode)
np <- map_dbl(data_sdm, function(x) nrow(x)/2)
ab_range <- c("US", "West", "West", "US", "US", "West", "West", "East",
              "East", "US", "East", "West", "US", "US", "US", "US", "US")

table1 <- data.frame(np = np[idx],
                     auc_hmsc_full = round(auc_hmsc_full[idx],2),
                     auc_brt_full = round(auc_brt_full[idx],2))

write.csv(table1, file = "table1a.csv")



table2a <- data.frame(spcode = spcode[idx],
                      prob_hmsc = round(beta_prob[[1]][idx], 3),
                      prob_brt = round(beta_prob_brt[[1]][idx], 3),
                      Rsq_hmsc = round(R_sq[[1]][idx], 2),
                      Rsq_brt = round(R_sq_brt[[1]][idx], 2))

write.csv(table2a, file = "table2a.csv")

table2b <- data.frame(spcode = spcode[idx],
                      prob_hmsc = round(beta_prob[[2]][idx], 3),
                      prob_brt = round(beta_prob_brt[[2]][idx], 3),
                      Rsq_hmsc = round(R_sq[[2]][idx], 2),
                      Rsq_brt = round(R_sq_brt[[2]][idx], 2))

write.csv(table2b, file = "table2b.csv")




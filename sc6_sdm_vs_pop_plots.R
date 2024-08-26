
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

theme_set(theme_bw())

results_brt <- readRDS("results/results_brt_ebird_200.rds")
results_hmsc_maps <- readRDS("results/results_hmsc_maps_200.rds")
results_brt_maps <- readRDS("results/results_brt_maps_200.rds")
results_cjspop <- readRDS("results/results_cjspop.rds")
results_dem <- readRDS("results/results_dem.rds")
results_pR <- readRDS("results/results_pR.rds")
auc_hmsc_full <- readRDS("results/results_auc_hmsc_200.rds")
auc_brt_full <- readRDS("results/results_auc_brt_200.rds")
data_sdm <- readRDS("data/data_sdm_ebird_pca_200.rds")
occ <- readRDS("data/data_sdm_ebird_raw_200.rds") 
occ <- map(occ, function(x) { 
  dplyr::select(x, y, Longitude, Latitude) %>%
    filter(y == 1) %>%
    dplyr::select(-y)
})
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

x_std <- purrr::map(results_hmsc_maps, function(x) (x - mean(x))/sd(x))

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
beta_chains[[1]] <- purrr::map(RvsPr, function(x) 
  MCMCchains(x, params = "beta"))
beta_chains[[2]] <- purrr::map(NvsPr, function(x) 
  MCMCchains(x, params = "beta"))
beta_chains[[3]] <- purrr::map(hNvsPr, function(x) 
  MCMCchains(x, params = "beta"))
beta_chains[[4]] <- purrr::map(pRvsPr, function(x) 
  MCMCchains(x, params = "beta"))

beta_prob <- purrr::map(beta_chains, function(x) 
  map_dbl(x, function(y) length(which(y>0))/length(y)))

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
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits = c(0,0.5)) +
  labs(y = bquote("R"^2), title = "GLM") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


# Repeat analysis with BRT results ----------------------------------------

x_std_brt <- purrr::map(results_brt_maps, function(x) (x - mean(x))/sd(x))

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
beta_chains_brt[[1]] <- purrr::map(RvsPr_brt, function(x) 
  MCMCchains(x, params = "beta"))
beta_chains_brt[[2]] <- purrr::map(NvsPr_brt, function(x) 
  MCMCchains(x, params = "beta"))
beta_chains_brt[[3]] <- purrr::map(hNvsPr_brt, function(x) 
  MCMCchains(x, params = "beta"))
beta_chains_brt[[4]] <- purrr::map(pRvsPr_brt, function(x) 
  MCMCchains(x, params = "beta"))

beta_prob_brt <- purrr::map(beta_chains_brt, function(x) 
  map_dbl(x, function(y) length(which(y>0))/length(y)))

df3 <- data.frame(beta = unlist(beta_prob_brt),
                  type = c(rep("R", length(beta_chains_brt[[1]])),
                           rep("N", length(beta_chains_brt[[2]])),
                           rep("highN", length(beta_chains_brt[[3]])),
                           rep("PR", length(beta_chains_brt[[4]]))))

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
R_sq_brt[[1]] <- map_dbl(RvsPr_brt, function(x) 
  MCMCsummary(x, params = "Rsq")[1,1])
R_sq_brt[[2]] <- map_dbl(NvsPr_brt, function(x) 
  MCMCsummary(x, params = "Rsq")[1,1])
R_sq_brt[[3]] <- map_dbl(hNvsPr_brt, function(x) 
  MCMCsummary(x, params = "Rsq")[1,1])
R_sq_brt[[4]] <- map_dbl(pRvsPr_brt, function(x) 
  MCMCsummary(x, params = "Rsq")[1,1])

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
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits = c(0,0.5)) +
  labs(y = bquote("R"^2), title = "BRT") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10))


g3 + g1 + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("figures/fig1.jpeg", width = 16.8, height = 12, units = "cm")

g4 + g2 + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("figures/fig2.jpeg", width = 16.8, height = 12, units = "cm")


# Comparison of Pocc with r<0 and r>0 -------------------------------------

# Select species with negative r populations
neg_r <- map_dbl(popR, function(x) sum(log(x)<0))
pos_r <- map_dbl(popR, function(x) sum(log(x)>0))

ss <- apply(cbind(neg_r, pos_r), 1, function(x) paste(x[2], x[1], sep=" - ")) 

popR2 <- popR[neg_r>=2]
pocc_hmsc <- results_hmsc_maps[neg_r>=2]
pocc_brt <- results_brt_maps[neg_r>=2]

r_group <- map(popR2, function(x) as.factor(as.numeric(log(x)>0)))

auc_brt <- map2_dbl(r_group, pocc_brt, 
                    function(x,y) auc(x, y, direction = "<"))
auc_hmsc <- map2_dbl(r_group, pocc_hmsc, 
                     function(x,y) auc(x, y, direction = "<"))

spcode2 <- spcode[neg_r>=2]
spcode3 <- map2(spcode2, map_dbl(popR2, length), function(x,y) rep(x, y)) %>%
  unlist()
spcode3 <- factor(spcode3, levels = spcode2[order(auc_brt, decreasing = T)])

box_data1 <- data.frame(y = unlist(pocc_hmsc),
                        x = spcode3,
                        group = factor(unlist(r_group),
                                       levels = c(1,0)))

box_data2 <- data.frame(y = unlist(pocc_brt),
                        x = spcode3,
                        group = factor(unlist(r_group),
                                       levels = c(1,0)))

g5 <- ggplot() +
  geom_boxplot(mapping = aes(y = y, x = x, fill = group), 
               data = box_data1, alpha = 0.7) +
  geom_text(mapping = aes(x = spcode2, y = 1.03, 
                          label = round(auc_hmsc, 2)), size = 3, 
            fontface = "bold") +
  geom_text(mapping = aes(x = spcode2, y = -0.05, 
                          label = ss[neg_r>=2]), size = 3) +
  scale_fill_manual(labels = c(bquote(bar("r")*">0"), bquote(bar("r")*"<0")), 
                    values = c("orange", "grey")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Suitability", title = "GLM")

g6 <- ggplot() +
  geom_boxplot(mapping = aes(y = y, x = x, fill = group), 
               data = box_data2, alpha = 0.7) +
  geom_text(mapping = aes(x = spcode2, y = 1.03, 
                          label = round(auc_brt,2)), size = 3,
            fontface = "bold") +
  geom_text(mapping = aes(x = spcode2, y = -0.05, 
                          label = ss[neg_r>=2]), size = 3) +
  scale_fill_manual(labels = c(bquote(bar("r")*">0"), bquote(bar("r")*"<0")), 
                    values = c("orange", "grey")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Suitability", title = "BRT")

g6/g5 +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag.position = c(0.05, 1))
ggsave("figures/fig3.jpeg", width = 20, height = 16.8, units = "cm")


# Maps --------------------------------------------------------------------

library(maps, pos = 30)
library(ggspatial)

data_midpoint <- readRDS("data/data_midpoint.rds")

usa <- st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))
laea <- st_crs("+proj=laea +lat_0=30 +lon_0=-95") # Lambert equal area
usa <- st_transform(usa, laea)

draw_rmap <- function(i, results_sdm, lgd.pos = "bottom",
                      xmax = -1000000,
                      ymax = -280000,
                      pad_x = 1,
                      pad_y = 0.7) {
  
  pop_no <- row.names(results_cjspop$chdata[[i]]$Nobs_ad)
  pop_coord <- filter(data_midpoint, pop %in% pop_no)
  pop_coord <- pop_coord[order(pop_coord$pop),]
  
  suit <- results_sdm[[i]] 
  
  occ_dat <- occ[[i]]
  occ_dat <- st_multipoint(as.matrix(occ_dat)) %>%
    st_sfc()
  st_crs(occ_dat) <- 4326
  occ_dat <- st_transform(occ_dat, laea)
  
  z <- foreach(h = 1:nrow(pop_coord)) %do% {
    st_point(as.matrix(pop_coord[h,3:4]))
  }
  z_sfc <- st_sfc(z)
  st_crs(z_sfc) <- 4326
  z_sfc <- st_transform(z_sfc, laea)
  
  sp_index <- which(names(r_group) == as.character(spcode[i]))
  
  word1 <- sp_eng[i]
  word2 <- as.character(spcode[i])
  word3 <- sp_lat[i]
  
  if (length(sp_index > 0)) {
    d <- st_sf(data.frame(group = r_group[[sp_index]], 
                          Suitability = suit, 
                          geom = z_sfc))
    
    ggplot() +
      geom_sf(data = usa, alpha = 0.9, fill = "snow1") +
      geom_sf(data = occ_dat, alpha = 0.5, shape = ".") +
      geom_sf(data = d, size = 2, alpha = 0.9, 
              aes(shape = group, color = Suitability)) +
      #scale_color_gradientn(colours = terrain.colors(10, rev = T),
                            #limits = c(0,1)) +
      scale_colour_gradient2_tableau(
        palette = "Orange-Blue Diverging", limits = c(0,1)) +
      scale_shape_manual(values = c(17, 19), 
                         labels = c("sink", "source"),
                         name = "MAPS\nLocations") +
      labs(title = bquote(.(word1)*~"("*.(word2)*","~italic(.(word3))*")")) +
      annotate("rect", xmin = -2450853.4, xmax = xmax, 
               ymin = -457753.3, ymax = ymax, fill = "white") +
      annotation_scale(mapping = aes(location = "bl"),
                       pad_x = unit(pad_x, "cm"),
                       pad_y = unit(pad_y, "cm")) +
      theme(legend.position = lgd.pos,
            panel.background = element_rect(fill = "lightblue")) +
      scale_y_continuous(limits = c(-457753.3, 2381225.5 )) +
      scale_x_continuous(limits = c(-2450853.4, 2186391.9))
  } else {
    d <- st_sf(data.frame(Suitability = suit, 
                          geom = z_sfc))
    
    ggplot() +
      geom_sf(data = usa, alpha = 0.9, fill = "snow1") +
      geom_sf(data = occ_dat, alpha = 0.5, shape = ".") +
      geom_sf(data = d, size = 2, alpha = 0.9, 
              aes(color = Suitability)) +
      #scale_color_gradientn(colours = terrain.colors(10, rev = T)) +
      scale_colour_gradient2_tableau(
        palette = "Orange-Blue Diverging", limits = c(0,1)) +
      labs(title = bquote(.(word1)*~"("*.(word2)*","~italic(.(word3))*")")) +
      annotate("rect", xmin = -2450853.4, xmax = -1000000, 
               ymin = -457753.3, ymax = -280000, fill = "white") +
      annotation_scale(mapping = aes(location = "bl"),
                       pad_x = unit(1, "cm"),
                       pad_y = unit(0.7, "cm")) +
      theme(legend.position = lgd.pos,
            panel.background = element_rect(fill = "lightblue")) +
      scale_y_continuous(limits = c(-457753.3, 2381225.5 )) +
      scale_x_continuous(limits = c(-2450853.4, 2186391.9))
  }
}

sp_lat <- c("Aimophila ruficeps", "Vireo huttoni","Melozone crissalis",
            "Dryobates villosus", "Certhia americana", "Contopus sordidulus",
            "Poecile gambeli", "Poecile carolinensis", "Empidonax virescens",
            "Melospiza lincolnii", "Baeolophus bicolor", "Chamaea fasciata",
            "Dryobates pubescens", "Vireo olivaceus", "Passerina cyanea",
            "Pipilo maculatus", "Icteria virens")

sp_eng <- c("Rufous-crowned Sparrow", "Hutton's Vireo", "California Towhee", 
            "Hairy Woodpecker", "Brown Creeper", "Western Wood-pewee", 
            "Mountain Chickadee", "Carolina Chickadee", "Acadian Flycatcher", 
            "Lincoln's Sparrow", "Tufted Titmouse", "Wrentit", 
            "Downy Woodpecker", "Red-eyed Vireo", "Indigo Bunting", 
            "Spotted Towhee", "Yellow-breasted Chat")

g_map1 <- draw_rmap(16, results_brt_maps, lgd.pos = "none",
                    xmax = -750000, ymax = -240000, pad_x = 0.6, pad_y = 0.39)
g_map2 <- draw_rmap(4, results_brt_maps, lgd.pos = "none",
                    xmax = -750000, ymax = -240000, pad_x = 0.6, pad_y = 0.39)

draw_rmap(4, results_brt_maps, lgd.pos = "right")
ggsave("figures/map_legend.pdf", width = 20, height = 20, units = "cm")

for (i in 1:17) {
  print(draw_rmap(i, results_brt_maps, lgd.pos = "bottom"))
  ggsave(paste("figures/sup_map", i, ".jpeg", sep = ""),
         width = 20, height = 20, units = "cm")
}


# r, N vs Pocc plots ------------------------------------------------------

plot_r <- function(i, gmap) {
  
  n_chains <- MCMCchains(RvsPr_brt[[i]], params = c("alpha", "beta"))
  y <- list()
  for (h in 1:length(x_std[[i]])) {
    
    y[[h]] <- n_chains[,1] + n_chains[,2]*x_std_brt[[i]][h]
    
  }
  
  y_avg <- map_dbl(y, mean)
  y_low <- map_dbl(y, function(x) quantile(x, 0.025))
  y_high <- map_dbl(y, function(x) quantile(x, 0.975))
  
  data1 <- data.frame(x = results_brt_maps[[i]],
                      y = log(popR[[i]]))
  
  data2 <- data.frame(x = results_brt_maps[[i]],
                      y = y_avg,
                      y_low = y_low,
                      y_high = y_high)
  
  ggplot(data = data1, aes(x = x, y = y)) +
    geom_point(aes(color = gmap$layers[[3]]$data$Suitability,
                   shape = gmap$layers[[3]]$data$group),
               size = 2) +
    geom_line(data = data2, aes(x = x, y = y), linewidth = 1.5,
              color = "darkgoldenrod4", linetype = 2) +
    geom_ribbon(data = data2, aes(ymin = y_low, ymax = y_high), alpha = 0.4, 
                fill = "grey") +
    scale_shape_manual(values = c(17, 19), 
                       labels = c("sink", "source")) +
    scale_colour_gradient2_tableau(
      palette = "Orange-Blue Diverging", limits = c(0,1)) +
    labs(x = "Suitability", 
         y = bquote("Intrinsic growth rate"*~"("*bar("r")*")")) +
    theme(legend.position = "none",
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 12))
}

gr1 <- plot_r(16, g_map1)
gr2 <- plot_r(4, g_map2)

(g_map1 | gr1) /
  (g_map2 | gr2) +
  plot_annotation(tag_levels = "a")

ggsave("figures/fig4.pdf", width = 25, height = 20, units = "cm")  

## Intrinsic growth
rg_brt <- list()
for (i in 1:17) {
  
  n_chains <- MCMCchains(RvsPr_brt[[i]], params = c("alpha", "beta"))
  y <- list()
  for (h in 1:length(x_std[[i]])) {
    
    y[[h]] <- n_chains[,1] + n_chains[,2]*x_std_brt[[i]][h]
    
  }
  
  y_avg <- map_dbl(y, mean)
  y_low <- map_dbl(y, function(x) quantile(x, 0.025))
  y_high <- map_dbl(y, function(x) quantile(x, 0.975))
  
  data1 <- data.frame(x = results_brt_maps[[i]],
                      y = log(popR[[i]]))
  
  data2 <- data.frame(x = results_brt_maps[[i]],
                      y = y_avg,
                      y_low = y_low,
                      y_high = y_high)
  
  auc_sp <- auc_brt[which(names(auc_brt) == spcode[i])]
  if (length(auc_sp) < 1) auc_sp <- NA
  
  title <- paste(sp_eng[i], " (", "BRT, ", "AUC = ", 
                 round(auc_sp, 2), ")", sep = "")
  
  rg_brt[[i]] <- ggplot(data = data1, aes(x = x, y = y)) +
    geom_point() +
    geom_line(data = data2, aes(x = x, y = y), color = "orange") +
    geom_ribbon(data = data2, aes(ymin = y_low, ymax = y_high), alpha = 0.2, 
                fill = "orange") +
    labs(x = "Suitability", y = bquote(bar("r")), title = title) +
    theme(plot.title = element_text(size = 6),
          axis.title = element_text(size = 8))
  
}

rg_brt <- rg_brt[order(spcode)]

## Abundance
ng_brt <- list()
for (i in 1:17) {
  
  n_chains <- MCMCchains(NvsPr_brt[[i]], params = c("alpha", "beta"))
  y <- list()
  for (h in 1:length(x_std_brt[[i]])) {
    
    y[[h]] <- n_chains[,1] + n_chains[,2]*x_std_brt[[i]][h]
    
  }
  
  y_avg <- map_dbl(y, mean)
  y_low <- map_dbl(y, function(x) quantile(x, 0.025))
  y_high <- map_dbl(y, function(x) quantile(x, 0.975))
  
  data1 <- data.frame(x = results_brt_maps[[i]],
                      y = log(N[[i]]$popN))
  
  data2 <- data.frame(x = results_brt_maps[[i]],
                      y = y_avg,
                      y_low = y_low,
                      y_high = y_high)
  
  auc_sp <- auc_brt[which(names(auc_brt) == spcode[i])]
  if (length(auc_sp ) < 1) auc_sp <- NA
  
  title <- paste(sp_eng[i], " (", "BRT, ", "AUC = ", 
                 round(auc_sp, 2), ")", sep = "")
  
  ng_brt[[i]] <- ggplot(data = data1, aes(x = x, y = y)) +
    geom_point() +
    geom_line(data = data2, aes(x = x, y = y), color = "orange") +
    geom_ribbon(data = data2, aes(ymin = y_low, ymax = y_high), alpha = 0.2, 
                fill = "orange") +
    labs(x = "Suitability", y = bquote(bar("N")), title = title) +
    theme(plot.title = element_text(size = 6),
          axis.title = element_text(size = 8))
  
}

ng_brt <- ng_brt[order(spcode)]

(ng_brt[[1]] + rg_brt[[1]]) / 
  (ng_brt[[2]] + rg_brt[[2]]) / 
  (ng_brt[[3]] + rg_brt[[3]])  +
  plot_annotation(tag_levels = 'a')
ggsave("figures/supp_fig_r1.jpeg", width = 6, height = 6, units = "in")

(ng_brt[[4]] + rg_brt[[4]]) / 
  (ng_brt[[5]] + rg_brt[[5]]) / 
  (ng_brt[[6]] + rg_brt[[6]])  +
  plot_annotation(tag_levels = 'a')
ggsave("figures/supp_fig_r2.jpeg", width = 6, height = 6, units = "in")

(ng_brt[[7]] + rg_brt[[7]]) / 
  (ng_brt[[8]] + rg_brt[[8]]) / 
  (ng_brt[[9]] + rg_brt[[9]]) +
  plot_annotation(tag_levels = 'a')
ggsave("figures/supp_fig_r3.jpeg", width = 6, height = 6, units = "in")

(ng_brt[[10]] + rg_brt[[10]]) / 
  (ng_brt[[11]] + rg_brt[[11]]) / 
  (ng_brt[[12]] + rg_brt[[12]]) +
  plot_annotation(tag_levels = 'a')
ggsave("figures/supp_fig_r4.jpeg", width = 6, height = 6, units = "in")

(ng_brt[[13]] + rg_brt[[13]]) / 
  (ng_brt[[14]] + rg_brt[[14]]) / 
  (ng_brt[[15]] + rg_brt[[15]]) +
  plot_annotation(tag_levels = 'a')
ggsave("figures/supp_fig_r5.jpeg", width = 6, height = 6, units = "in")

(ng_brt[[16]] + rg_brt[[16]]) / (ng_brt[[17]] + rg_brt[[17]]) +
  plot_annotation(tag_levels = 'a')
ggsave("figures/supp_fig_r6.jpeg", width = 6, height = 6, units = "in")


# Tables ------------------------------------------------------------------

idx <- order(spcode)
np <- map_dbl(data_sdm, function(x) nrow(x)/2)

table1 <- data.frame(np = np[idx],
                     auc_hmsc_full = round(auc_hmsc_full[idx],2),
                     auc_brt_full = round(auc_brt_full[idx],2))

write.csv(table1, file = "tables/table1a.csv")

table2a <- data.frame(spcode = spcode[idx],
                      prob_hmsc = round(beta_prob[[1]][idx], 3),
                      prob_brt = round(beta_prob_brt[[1]][idx], 3),
                      Rsq_hmsc = round(R_sq[[1]][idx], 2),
                      Rsq_brt = round(R_sq_brt[[1]][idx], 2))

write.csv(table2a, file = "tables/table2a.csv")

table2b <- data.frame(spcode = spcode[idx],
                      prob_hmsc = round(beta_prob[[2]][idx], 3),
                      prob_brt = round(beta_prob_brt[[2]][idx], 3),
                      Rsq_hmsc = round(R_sq[[2]][idx], 2),
                      Rsq_brt = round(R_sq_brt[[2]][idx], 2))

write.csv(table2b, file = "tables/table2b.csv")


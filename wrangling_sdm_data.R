# Prepare occurence and climate PCA for SDMs

rm(list = ls())

library(tidyverse)
library(raster)
library(factoextra)


# Load occ data -----------------------------------------------------------

spname <- c("Aimophila ruficeps", "Vireo huttoni","Melozone crissalis",
            "Dryobates villosus", "Certhia americana", "Contopus sordidulus",
            "Poecile gambeli", "Poecile carolinensis", "Empidonax virescens",
            "Melospiza lincolnii", "Baeolophus bicolor", "Chamaea fasciata",
            "Dryobates pubescens", "Vireo olivaceus", "Passerina cyanea",
            "Pipilo maculatus", "Icteria virens")

spcode <- readRDS("results_cjspop.rds")$spcode

# Occurences from Gbif
occ_temp <- readRDS("data_occ.rds")
occ <- list()

# Remove duplicates between 1992-2008
for (i in 1:length(spname)) {
  occ[[i]] <- filter(occ_temp, grepl(spname[i], species)) %>%
    filter(lat >=20 & lat <=50 & long >=-140 & long<=-50) %>%
    mutate_at(.vars = "long", function(x) x+360) %>%
    dplyr::select(long,lat) %>%
    as.matrix(.)
}
names(occ) <- spcode

# Load rasters ------------------------------------------------------------

rasters <- paste("weather/average", list.files("weather/average"), sep = "/")
r_avg <- stack(rasters)

# Create SDM data ---------------------------------------------------------

# Filter presences
occ1 <- list()
occ2 <- list()
for (i in 1:length(occ)) {
  
  z <- raster::extract(r_avg, occ[[i]], cellnumbers = T) %>%
    as.data.frame() %>% 
    group_by(cells) %>% summarise(n = n())
  
  # Remove duplicates
  cells_index <- distinct(z, cells)$cells
  occ1[[i]] <- occ[[i]][z$cells %in% cells_index,]
  
  # Remove cells with only 1 observation
  occ2[[i]] <- occ[[i]][which(z$n>1),]
  
}


# Extract weather data and generate pseudo-absences
prep.pa <- function(occ,r) {
  
  p <- raster::extract(r, occ, cellnumbers = T) %>%
    cbind(occ, .)
  
  index <- which(is.na(p), arr.ind = T)[,1]
  if (length(index>0)) p <- p[-index,]
  
  # NA to presence points
  r2 <- r
  for (k in 1:length(r@layers)) r2[[k]][p[,"cells"]] <- NA
  
  # Sample pseudo-absences
  absence <- rasterToPoints(r2)
  if (nrow(p)<nrow(absence)) {
    
    index2 <- sample(1:nrow(absence), size = nrow(p))
    absence <- absence[index2,]
    
  }
  colnames(absence)[1:2] <- c("long", "lat")
  
  index3 <- which(is.na(absence), arr.ind = T)[,1]
  if (length(index3>0)) absence <- absence[-index3,]
  
  p <- p[,-3]
  pa <- rbind(p, absence) %>%
    as.data.frame() %>%
    mutate(y = c(rep(1, nrow(p)), rep(0, nrow(absence))))
  
  return(pa)
}

pb <- txtProgressBar(0, length(occ), style = 3)
data_sdm_raw1 <- list()
data_sdm_raw2 <- list()

for (i in 1:length(occ)) {
    
  data_sdm_raw1[[i]] <- prep.pa(occ = occ1[[i]], r = r_avg)
  data_sdm_raw2[[i]] <- prep.pa(occ = occ2[[i]], r = r_avg)
  
  setTxtProgressBar(pb,i)
  
}

names(data_sdm_raw1) <- spcode
names(data_sdm_raw2) <- spcode
saveRDS(data_sdm_raw1, file = "data_sdm_raw1.rds")
saveRDS(data_sdm_raw2, file = "data_sdm_raw2.rds")

# PCA on raw sdm data
apply.pca <- function(data_sdm, spcode) {
  
  res_pca <- map(data_sdm, function(x) prcomp(x[,3:12], scale. = TRUE))
  names(res_pca) <- spcode
  
  # Select dimensions that explain 80% of the variance
  eig <- map(res_pca, function(x) get_eigenvalue(x))
  dim_num <- map(eig, function(x) which(x$cumulative.variance.percent>80)[1])
  
  # PCA data for occurences and pseduo-absence
  y <- map(data_sdm, function(x) x$y)
  data_sdm_pca <- map2(res_pca, dim_num, function(x,z) get_pca_ind(x)$coord[,1:z]) %>%
    map2(., y, function(x,z) cbind(x,z) %>% as.data.frame())
  
  names(data_sdm_pca) <- spcode
  
  res <- list(res_pca = res_pca,
              data = data_sdm_pca)
  
  return(res)
  
}

data_sdm_pca1 <- apply.pca(data_sdm_raw1, spcode)
data_sdm_pca2 <- apply.pca(data_sdm_raw2, spcode)

saveRDS(data_sdm_pca1, file = "data_sdm_pca1.rds")
saveRDS(data_sdm_pca2, file = "data_sdm_pca2.rds")


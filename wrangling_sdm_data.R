# Prepare occurence and climate PCA for SDMs

library(tidyverse)
library(raster)
library(factoextra)
library(doSNOW)
library(foreach)
library(sf)

# Load occ data -----------------------------------------------------------

spname <- c("Aimophila ruficeps", "Vireo huttoni","Melozone crissalis",
            "Dryobates villosus", "Certhia americana", "Contopus sordidulus",
            "Poecile gambeli", "Poecile carolinensis", "Empidonax virescens",
            "Melospiza lincolnii", "Baeolophus bicolor", "Chamaea fasciata",
            "Dryobates pubescens", "Vireo olivaceus", "Passerina cyanea",
            "Pipilo maculatus", "Icteria virens")

spcode <- readRDS("results_cjspop.rds")$spcode

chdata <- readRDS("results_cjspop.rds")$chdata

# Occurences from Gbif
occ_temp <- readRDS("data_occ.rds")
occ <- list()

# Remove duplicates between 1992-2008
for (i in 1:length(spname)) {
  occ[[i]] <- filter(occ_temp, grepl(spname[i], species)) %>%
    filter(lat >=20 & lat <=50 & long >=-140 & long<=-50) %>%
    dplyr::select(long,lat) %>%
    as.matrix(.)
}
names(occ) <- spcode


# Load rasters and polygons -----------------------------------------------

rasters <- paste("weather/average", list.files("weather/average"), sep = "/")
r_avg <- stack(rasters)
r_avg <- rotate(r_avg)

physio <- st_read("physio_shp/physio.shp")

idx <- list()
for (i in 1:25) {
  
  if (is.na(unique(physio$PROVINCE)[i])) {
    
    idx[[i]] <- NA
    
  } else {
    
    idx[[i]] <- which(physio$PROVINCE %in% unique(physio$PROVINCE)[i])
    
  }
  
}
idx <- idx[-18]

## Range for generating pseudo absences
west <- physio[unlist(idx[c(3,4,5,7,9,13,15,16,19,24)]),8]
east <- physio[unlist(idx[c(1,2,5,6,8,9,10,11,12,13,14,16,17,18,20:23)]),8]
US <- physio

ab_range <- list(US, west, west, US, US, west, west, east,
                 east, US, east, west, US, US, US, US, US)

names(ab_range) <- spcode

# Create SDM data ---------------------------------------------------------

# Filter presences
occ_filter <- list()
for (i in 1:length(occ)) {
  
  z <- raster::extract(r_avg, occ[[i]], cellnumbers = T) %>%
    as.data.frame() %>% 
    group_by(cells) %>% summarise(n = n())
  
  # Remove cells with only 1 observation
  occ_filter[[i]] <- occ[[i]][which(z$n>1),]
  
}

names(occ_filter) <- spcode

saveRDS(occ_filter, file = "data_occ_filter.rds")

# Extract weather data and generate pseudo-absences
prep.pa <- function(occ, r, ab_no, ab_range) {
  
  p <- raster::extract(r, occ, cellnumbers = T) %>%
    cbind(occ, .)
  
  index <- which(is.na(p), arr.ind = T)[,1]
  if (length(index>0)) p <- p[-index,]
  
  r2 <- mask(r, ab_range)
  
  # NA to presence points
  for (k in 1:length(r@layers)) r2[[k]][p[,"cells"]] <- NA
  
  # Sample pseudo-absences
  absence <- rasterToPoints(r2)
  
  if (ab_no == "equal") {
    
    index2 <- sample(1:nrow(absence), size = nrow(p))
    absence <- absence[index2,]
    
  } else {
    
    index2 <- sample(1:nrow(absence), size = ab_no)
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

data_sdm_raw <- list()
pb <- txtProgressBar(0,17,style = 3)
for (i in 1:17) {
  
  data_sdm_raw[[i]] <- prep.pa(occ = occ_filter[[i]], 
                               r = r_avg, 
                               ab_range = ab_range[[i]], 
                               ab_no = "equal")
  setTxtProgressBar(pb,i)
  
}

names(data_sdm_raw) <- spcode
saveRDS(data_sdm_raw, file = "data_sdm_raw.rds")

# PCA on raw sdm data
apply.pca <- function(data_sdm, res_pca) {
  
  # Select dimensions that explain at least 80% of the variance
  eig <- get_eigenvalue(res_pca)
  dim_num <- which(eig$cumulative.variance.percent>80)[1]
  
  # PCA data for occurences and pseduo-absence
  y <- data_sdm$y
  data_sdm_pca <- predict(res_pca, data_sdm)[,1:dim_num] %>%
    cbind(.,y) %>%
    as.data.frame()
  
  return(data_sdm_pca)
  
}

data_sdm_pca <- foreach (i=1:17) %do% {

  res <- apply.pca(data_sdm_raw[[i]], chdata[[i]]$res_pca)
  
  return(res)
  
}
names(data_sdm_pca) <- spcode

saveRDS(data_sdm_pca, file = "data_sdm_pca.rds")


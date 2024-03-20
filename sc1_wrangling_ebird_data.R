
# Prepare occurrence data from eBird and climate PCA for SDMs

library(tidyverse)
library(raster)
library(factoextra)
library(doSNOW)
library(foreach)
library(sf)
library(sp)
library(auk)
library(spThin)
library(dismo)

spcode <- readRDS("results/results_cjspop.rds")$spcode
chdata <- readRDS("results/results_cjspop.rds")$chdata


# read in and filter eBird data -------------------------------------------

extract_ebird <- function(spcode) {
  directory <- "/Volumes/T7 Shield/ebird/"
  ebd_folder <- paste(spcode, "/", sep = "")
  ebd_file <- paste(spcode, "ebd.txt", sep = "_")
  
  ebird_folder <- "ebird/"
  
  f_in <- paste(directory, ebd_folder, ebd_file, sep = "")
  f_out <- paste(ebird_folder, spcode, "_filter.txt", sep = "")
  
  ebird_data <- f_in %>% 
    # 1. reference file
    auk_ebd() %>% 
    # 2. define filters
    auk_country("US") %>%
    auk_date(date = c("*-05-01", "*-08-31")) %>%
    auk_bbox(bbox = c(-140, 20, -50, 50)) %>%
    auk_time(start_time = c("04:15", "12:00")) %>%
    auk_duration(duration = c(0, 180)) %>%
    auk_distance(c(0, 3.1)) %>%
    auk_complete() %>%
    # 3. run filtering
    auk_filter(file = f_out, overwrite = T) %>% 
    # 4. read text file into r data frame
    read_ebd()
  
  ebird_data$year <- year(ebird_data$observation_date)
  
  ebird_data <- filter(ebird_data, year >= 1992 & year <= 2008) %>%
    dplyr::select(longitude, latitude)
  
  return(ebird_data)
   
}

pb <- txtProgressBar(0, 17, style = 3)

occ_ebird <- foreach(i = 1:17) %do% {
  
  setTxtProgressBar(pb, i)
  
  extract_ebird(spcode = spcode[i])
}

## Thin occurrence records by 10 km
thin.occ <- function(occ) {
  occ$SPEC <- "x"
  z <- thin(occ, lat.col = "latitude", long.col = "longitude", reps = 10,
            thin.par = 10, locs.thinned.list.return = T, write.files = F)
  data_length <- map_dbl(z, nrow)
  idx <- which(data_length == max(data_length))
  idx2 <- sample(idx, 1)
  z[[idx2]]
}

thin.occ.long <- function(occ) {
  half_point <- round(nrow(occ)/2)
  occ1 <- occ[1:half_point,]
  occ2 <- occ[(half_point + 1):nrow(occ),]
  
  z1 <- thin.occ(occ1)
  z2 <- thin.occ(occ2)
  
  z3 <- rbind(z1, z2)
  colnames(z3) <- c("longitude", "latitude")
  
  thin.occ(z3)
}

occ_thin <- foreach(i = 1:17) %do% {
  if (nrow(occ_ebird[[i]]) > 20000) {
    thin.occ.long(occ_ebird[[i]])
  } else {
    thin.occ(occ_ebird[[i]])
  }
} 


# Load rasters and polygons -----------------------------------------------

rasters <- paste("weather/average", list.files("weather/average"), sep = "/")
r_avg <- stack(rasters)
r_avg <- rotate(r_avg)


# Create SDM data ---------------------------------------------------------

# Extract weather data and generate pseudo-absences
prep.pa <- function(occ, r, ab_no, buf_width) {
  
  p <- raster::extract(r, occ) %>%
    cbind(occ, .)
  
  # remove cells with NAs
  index <- which(is.na(p), arr.ind = T)[,1]
  if (length(index>0)) p <- p[-unique(index),]
  
  # Buffer to generate pseudoabsences
  buf1 <- raster::buffer(SpatialPoints(occ), width = buf_width)
  buf2 <- raster::buffer(SpatialPoints(occ), width = 10000)
  r2 <- mask(r, buf1)
  r3 <- mask(r2, buf2, inverse = T)
  
  absence <- randomPoints(r3, nrow(p))
  colnames(absence)[1:2] <- c("Longitude", "Latitude")
  a <- raster::extract(r, absence) %>%
    cbind(absence, .)
  
  # remove cells with NAs
  index <- which(is.na(a), arr.ind = T)[,1]
  if (length(index>0)) a <- a[-unique(index),]
  
  # Combine Presences and Absences 
  pa <- rbind(p, a) %>%
    as.data.frame() %>%
    mutate(y = c(rep(1, nrow(p)), rep(0, nrow(a))))
  
  return(pa)
}

data_sdm_raw_500 <- foreach(i = 1:17) %do% {
  prep.pa(occ = occ_thin[[i]], 
          r = r_avg, 
          buf_width = 500000)
}

data_sdm_raw_200 <- foreach(i = 1:17) %do% {
  prep.pa(occ = occ_thin[[i]], 
          r = r_avg, 
          buf_width = 200000)
}

data_sdm_raw_1000 <- foreach(i = 1:17) %do% {
  prep.pa(occ = occ_thin[[i]], 
          r = r_avg, 
          buf_width = 1000000)
}

names(data_sdm_raw_500) <- spcode
names(data_sdm_raw_200) <- spcode
names(data_sdm_raw_1000) <- spcode

saveRDS(data_sdm_raw_500, file = "data_sdm_ebird_raw_500.rds")
saveRDS(data_sdm_raw_200, file = "data_sdm_ebird_raw_200.rds")
saveRDS(data_sdm_raw_1000, file = "data_sdm_ebird_raw_1000.rds")

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

data_sdm_pca_500 <- foreach (i=1:17) %do% {
  
  res <- apply.pca(data_sdm_raw_500[[i]], chdata[[i]]$res_pca)
  
  return(res)
  
}
names(data_sdm_pca_500) <- spcode

data_sdm_pca_200 <- foreach (i=1:17) %do% {
  
  res <- apply.pca(data_sdm_raw_200[[i]], chdata[[i]]$res_pca)
  
  return(res)
  
}
names(data_sdm_pca_200) <- spcode

data_sdm_pca_1000 <- foreach (i=1:17) %do% {
  
  res <- apply.pca(data_sdm_raw_1000[[i]], chdata[[i]]$res_pca)
  
  return(res)
  
}
names(data_sdm_pca_1000) <- spcode

saveRDS(data_sdm_pca_500, file = "data/data_sdm_ebird_pca_500.rds")
saveRDS(data_sdm_pca_200, file = "data/data_sdm_ebird_pca_200.rds")
saveRDS(data_sdm_pca_1000, file = "data/data_sdm_ebird_pca_1000.rds")


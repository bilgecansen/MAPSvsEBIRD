
rm(list = ls())

library(tidyverse)
library(raster)


# Average of 17 years between 1992-2008  ----------------------------------

years <- 1992:2008
variables <- list.files("weather/single_year/1992")[1:10] %>%
  gsub(".asc", "", .)

if (!"average" %in% list.files("weather")) dir.create("weather/average")

pb <- txtProgressBar(min = 0, max = length(variables), style = 3)
for (k in 1:length(variables)) {

  r <- list()
  
  for (i in 1:length(years)) {
    
    folder <- paste("weather/single_year", years[i], sep = "/")
    r[[i]] <- raster(paste(folder, list.files(folder)[k], sep = "/"))

  }#i
    
  rs <- stack(r) %>% mean(.)
  names(r) <- paste(variables[k], years[i], sep = "_")
    
  filename <- paste(variables[k], "asc", sep = ".") %>%
    paste("weather/average", ., sep = "/")
    
  writeRaster(rs, filename = filename, overwrite = T)
  
  setTxtProgressBar(pb, k)
  
}#k

# Save as data frame
files <- paste0("weather/average/", list.files("weather/average"))
r <- map(files, raster) %>%
  do.call(raster::stack, .)
names(r) <- gsub(".asc", "", list.files("weather/average"))
df <- rasterToPoints(r)
saveRDS(df, file = "data_maurer.rds")


# Every year and variable in an array -------------------------------------

a <- list()
for (i in 1:length(years)) {
  
  folder <- paste("weather/single_year", years[i], sep = "/")
  
  r <- paste(folder, list.files(folder)[1:10], sep = "/") %>%
    stack(.)
  
  names(r) <- gsub(".asc", "", list.files(folder))[1:10]
  
  a[[i]] <- rasterToPoints(r)[,3:12]
  
  gc()
}

a <- lapply(a, function(x) array(x, dim = c(1,53097,10)))
a2 <- do.call(abind, list(a, along = 1))
dimnames(a2) <- list(years = as.character(years), NULL, variables = names(r))

saveRDS(a2, file = "maurer_array.rds")





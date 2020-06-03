
# Packages and data -------------------------------------------------------

library(tidyverse)
library(foreach)
library(Hmsc)

results_hmsc <- readRDS("results_hmsc.rds")
chdata <- readRDS("results_cjspop.rds")$chdata
results_pca <- readRDS("data_sdm_pca2.rds")$res_pca
data_sdm_pca <- readRDS("data_sdm_pca2.rds")$data
data_maurer <- readRDS("data_maurer.rds")


# Check model performance -------------------------------------------------

pb <- txtProgressBar(0, length(results_hmsc), style = 3)

cv <- list()
pred <- list()
model_fit <- list()
for (i in 1:length(results_hmsc)) {
  
  cv[[i]] <- createPartition(results_hmsc[[i]], nfolds = 10)
  pred[[i]] <- computePredictedValues(results_hmsc[[i]],
                                      partition = cv[[i]],
                                      nParallel = 4,
                                      expected = F)
  model_fit[[i]] <- evaluateModelFit(results_hmsc[[i]], pred[[i]])
  
  setTxtProgressBar(pb, i)
}

# AUC values
auc_hmsc <- map_dbl(model_fit, function(x) x$AUC) %>% print()
median(auc_hmsc)
saveRDS(auc_hmsc, file = "results_auc.rds")


# Predict occurence prob of MAPS locations --------------------------------

# Average weather values
avg_weather <- map(chdata, function(x) map_dfc(x$weather[1:10], function(y) apply(y, 1, mean)))
#for (i in 1:length(avg_weather)) names(avg_weather[[i]]) <- paste("layer", 1:10, sep = ".")

# convert raw weather data to pca variables
maps_pca <- map2(results_pca, avg_weather, predict) %>%
  map2(., data_sdm_pca, function(x,y) x[, 1:(ncol(y)-1)])

# Predict occurence prob with hmsc
maps_predict <- list()

for (i in 1:length(results_hmsc)) {
  
  maps_predict[[i]] <- predict(results_hmsc[[i]], XData = as.data.frame(maps_pca[[i]]), expected = T) %>%
    do.call(cbind, .) %>%
    apply(., 1, mean)
  
}

saveRDS(maps_predict, "results_maps_hmsc.rds")


# Predict occurence prob of whole US --------------------------------------

#colnames(data_maurer)[3:12] <- paste("layer", 1:10, sep = ".")


# convert raw weather data to pca variables
us_pca <- map(results_pca, function(x) predict(x, data_maurer[,3:12])) %>%
  map2(., data_sdm_pca, function(x,y) x[, 1:(ncol(y)-1)])

index <- map(us_pca, function(x) which(is.na(x), arr.ind = T)[,1])
us_pca2 <- map2(us_pca, index, function(x,y) x[-y,] )
coordinates <- data_maurer[-index[[1]], 1:2]

us_predict <- list()
pb <- txtProgressBar(min = 0, max = length(results_hmsc), style = 3)

for (i in 1:length(results_hmsc)) {
  
  us_predict[[i]] <- predict(results_hmsc[[i]], XData = as.data.frame(us_pca2[[i]]), expected = T) %>%
    do.call(cbind, .) %>%
    apply(., 1, mean)
  
  setTxtProgressBar(pb, i)
  
}

us_predict2 <- map(us_predict, function(x) cbind(coordinates, x))

saveRDS(us_predict2, "results_us_hmsc.rds")



# Packages and data -------------------------------------------------------

library(tidyverse)
library(foreach)
library(Hmsc)
library(dismo)
library(gbm)
library(tictoc)

results_hmsc <- readRDS("results_hmsc.rds")
results_brt <- readRDS("results_brt.rds")
chdata <- readRDS("results_cjspop.rds")$chdata
results_pca <- map(chdata, function(x) x$res_pca)
data_sdm_pca <- readRDS("data_sdm_pca.rds")


# Predict occurence prob of MAPS locations --------------------------------

# Average weather values
avg_weather <- map(chdata, function(x) map_dfc(x$weather[1:10], function(y) apply(y, 1, mean)))

# convert raw weather data to pca variables
maps_pca <- foreach(i=1:length(results_pca)) %do% {
    
  res <- predict(results_pca[[i]], avg_weather[[i]])
  res <- res[,1:(ncol(data_sdm_pca[[i]])-1)]
  
  return(res)
    
}
  
# Predict occurence prob with hmsc
maps_predict_hmsc <- foreach (i = 1:length(results_hmsc)) %do% {
    
    res <- predict(results_hmsc[[i]], 
                   XData = as.data.frame(maps_pca[[i]]), 
                   expected = T) %>% 
      do.call(cbind, .) %>%
      apply(., 1, mean)
    
    return(res)
    
}

names(maps_predict_hmsc) <- names(data_sdm_pca)

saveRDS(maps_predict_hmsc, "results_hmsc_maps.rds")

# Predict occurence prob with brt
maps_predict_brt <- foreach(i=1:17) %do% {
  
  colnames(maps_pca[[i]]) <- str_replace(colnames(maps_pca[[i]]), "Dim.", "PC")
  
  predict.gbm(results_brt[[i]], 
              as.data.frame(maps_pca[[i]]), 
              n.trees = results_brt[[i]]$n.trees,
              type = "response")
  
}  

names(maps_predict_brt) <- names(data_sdm_pca)

saveRDS(maps_predict_brt, file = "results_brt_maps.rds")


# Check model performance -------------------------------------------------

tic()
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

toc()

# AUC values
auc_hmsc_full <- map_dbl(model_fit, function(x) x$AUC) %>% print()
saveRDS(auc_hmsc_full, file = "results_auc_hmsc.rds")

auc_brt_full <- map_dbl(results_brt, function(x) x$cv.statistics$discrimination.mean)
saveRDS(auc_brt_full, file = "results_auc_brt.rds")


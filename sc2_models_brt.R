
# Boosted Regression Trees on selected species

library(dismo)
library(gbm)
library(doSNOW)
library(foreach)
library(tictoc)
library(tidyverse)

data_sdm_500 <- readRDS("data/data_sdm_ebird_pca_500.rds")
data_sdm_200 <- readRDS("data/data_sdm_ebird_pca_200.rds")
data_sdm_1000 <- readRDS("data/data_sdm_ebird_pca_1000.rds")
chdata <- readRDS("results/results_cjspop.rds")$chdata


# Boosted Regression Trees ------------------------------------------------

run_brt <- function(data_sdm) {
  cl <- makeCluster(6, types = "SOCK")
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = 17, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  tic()
  
  results <- foreach(i=1:17, .packages = "dismo", .options.snow = opts) %dopar% 
    {
      
      res <- gbm.step(data = data_sdm[[i]], 
                      gbm.x = 1:(ncol(data_sdm[[i]])-1), 
                      gbm.y = ncol(data_sdm[[i]]), 
                      tree.complexity = 2, 
                      family = "bernoulli", 
                      learning.rate = 0.01, 
                      bag.fraction = 0.75, 
                      step.size = 250, 
                      n.trees = 1000,
                      max.trees = 15000)
      
      if(!is.null(res)) {
        
        return(res)
      
      } else {
        
        res <- gbm.step(data = data_sdm[[i]], 
                        gbm.x = 1:(ncol(data_sdm[[i]])-1), 
                        gbm.y = ncol(data_sdm[[i]]), 
                        tree.complexity = 2, 
                        family = "bernoulli", 
                        learning.rate = 0.001, 
                        bag.fraction = 0.75, 
                        step.size = 250, 
                        n.trees = 1000,
                        max.trees = 15000)
        
        return(res)
        
      }
    
    }
  
  toc()
  
  stopCluster(cl)
  
  names(results) <- names(data_sdm_500)
  
  return(results)
}

results_500 <- run_brt(data_sdm_500)
results_200 <- run_brt(data_sdm_200)
results_1000 <- run_brt(data_sdm_1000)

saveRDS(results_500, file = "results/results_brt_ebird_500.rds")
saveRDS(results_200, file = "results/results_brt_ebird_200.rds")
saveRDS(results_1000, file = "results/results_brt_ebird_1000.rds")

  
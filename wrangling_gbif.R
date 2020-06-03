
rm(list = ls())
gc()

library(rgbif)
library(tidyverse)

spname <- c("Aimophila ruficeps", "Vireo huttoni","Melozone crissalis",
            "Dryobates villosus", "Certhia americana", "Contopus sordidulus",
            "Poecile gambeli", "Poecile carolinensis", "Empidonax virescens",
            "Melospiza lincolnii", "Baeolophus bicolor", "Chamaea fasciata",
            "Dryobates pubescens", "Vireo olivaceus", "Passerina cyanea",
            "Pipilo maculatus", "Icteria virens")

taxon_key <- map_chr(spname, function(x) name_suggest(q = x, rank = "species")$key)

# Download occ for 1992-2008 (for use with Maurer Historical) -------------

occ_req <- occ_download(
  paste("taxonKey = ", paste(taxon_key, collapse = ","), sep = ""),
  "country = US",
  "hasCoordinate = true",
  "year >= 1992",
  "year <= 2008",
  "month >= 5",
  "month <= 8",
  user = "bilgecan",
  pwd = "Mercimek1109",
  email = "bilgecan.sen@gmail.com"
)

# Check download status
occ_download_meta(occ_req)

# Import occ data
occ_dat <- occ_download_get(occ_req, overwrite = T) %>%
  occ_download_import(., fill=FALSE, quote = "") %>%
  filter(basisOfRecord == "HUMAN_OBSERVATION" & 
         (coordinateUncertaintyInMeters < 10000 | is.na(coordinateUncertaintyInMeters))  &
         hasGeospatialIssues == F &
         datasetName == "EOD - eBird Observation Dataset") %>%
  dplyr::select(species = scientificName,
         long = decimalLongitude,
         lat = decimalLatitude,
         year = year)

saveRDS(occ_dat, file = "data_occ.rds")


# MAPSvsEBIRD

Code for the manuscript titled "The potential for species distribution models to distinguish source populations from sinks ". 

**sc1_wrangling_ebird_data.R** prepares raw ebird and raster data for SDMs. This is for code sharing only, because these raw files are too large to share. The following scripts can be run as we share the data products of this code under the data folder.

**sc2_models_brt.R** runs the BRT SDMs and saves the output in the results folder.

**sc3_models_hmsc.R** runs the HMSC SDMs and saves the output in the results folder.

**sc4_summarize_sdm.R** prepares the SDM results for comparions with demographic data from MAPS locations.

**sc5_summarize_pR.R** calculates probability that a population is a sink (r<0) or source (r>0).

**sc6_sdm_vs_pop_plots.R** runs all analysis related to SDM vs demography comparisons and produces all the plots in the main manuscript. 

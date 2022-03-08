# MAPSvsGBIF

Code for the manuscript titled "Species distribution models can distinguish sources from sinks but fail to predict species demography". 

**wrangling_sdm_data.R** prepares data for BRT and HMSC models. We can't share the original data that goes into this script but we are sharing the output of it which are used to build these SDMs.

**models_brt.R** and **models_hmsc.R** runs the SDMs and **summarize_sdm.R** prepares the SDM results for further analysis.

**summarize_pR.R** calculates probability that a population is a sink (r<1) or source (r>1).

**sdm_vs_pop_plots** produces all the plots in the main manuscript.

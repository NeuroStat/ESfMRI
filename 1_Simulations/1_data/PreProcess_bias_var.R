####################
#### TITLE:     Preprocessing of the bias-variance Monte-Carlo simulations.
####            
#### Contents:
####
#### Source Files:
#### First Modified: 05/04/2018
#### Notes:
#################



##
###############
### Notes
###############
##

# This script takes the raw data of all Monte-Carlo simulations.
# Loads them into working memory, then saves one data frame.


##
###############
### Data
###############
##


# Provide data location
DataLoc_list <- list("TwoStageModels" = 
              "/Volumes/Elements/ReviewESfMRI/Simulations/TwoStageModels")
dataWD <- "TwoStageModels"
DataLoc <- DataLoc_list[[dataWD]]

# Empty data frame with results
VarHedgeRes <- tibble(sim = integer(),
                      param = factor(levels = c('beta1', 'S2G', 'hedge', 
                                  'varhedge', 'varhedge_radua',
                                  'beta1TR', 'S2GTR', 'hedgeTR', 'varhedgeTR')),
                      value = integer(),
                      nscan = numeric(),
                      nsub = numeric(),
                      trueValue = integer())

# Number of simulations
nsim <- 1000

# For loop over the simulations
for(i in 1:nsim){
  # Read in data object and bind to data frame
  VarHedgeRes <- bind_rows(VarHedgeRes,
                  readRDS(paste0(DataLoc, '/VarHedgeRes_BSvar_', i, '.rda')))
}

# Save data frame
saveRDS(object = VarHedgeRes, file = 'bias_variance_ES_fMRI.rda')


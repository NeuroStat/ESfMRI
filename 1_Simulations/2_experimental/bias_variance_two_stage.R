####################
#### TITLE:     Two-stage standardized effect approach 
#### Contents:
####
#### Source Files:
#### First Modified: 04/04/2018
#### Notes:
#################



##
###############
### Notes
###############
##


######## STOPPED WORKING ON THIS SCRIPT ########

# As we noticed in the bias_variance.R file, the estimation of the variance of
# a standardized effect size at second level fMRI is biased. 
# In this script, we consider a two-stage approach on the standardized effect sizes.
# That is: first calculate the ES for each subject. Then calculate a weighted average
# at group level. Its variance is then the inverse of the sum of the weights.


##
###############
### Preparation
###############
##


# Take argument from master file
input <- commandArgs(TRUE)
# K'th simulation
K <- try(as.numeric(as.character(input)[1]),silent=TRUE)

# Which machine
MACHINE <- try(as.character(input)[2],silent=TRUE)
# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
  K <- 1
}

# Set starting seed
starting.seed <- pi*K
set.seed(starting.seed)

# Set WD: this is location where results are written
if(MACHINE=='HPC'){
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/BiasVarTwoS/Results'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/Simulation/Results/BiasVar'
}

# Load in libraries
library(lattice)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(reshape2)
library(RColorBrewer)
library(mvmeta)
library(metafor)
library(devtools)
library(neuRosim)
library(NeuRRoStat)

# Data frame with results
# VarHedgeRes <- tibble(sim = integer(),
#                       param = factor(levels = c('beta1', 'S2G', 'hedge', 'varhedge', 
#                               'varhedge_radua',
#                               'beta1TR', 'S2GTR', 'hedgeTR', 'varhedgeTR')),
#                       value = integer(),
#                       nscan = numeric(),
#                       nsub = numeric(),
#                       trueValue = integer())
# 

##
###############
### Functions
###############
##

##
###############
### Simulation parameters
###############
##

# Vector with number of scans
Vnscan <- seq(100, 500, by = 10)

# Same with number of subjects
Vnsub <- seq(20, 100, by = 10)

# Combinations
combPar <- expand.grid('nscan' = Vnscan, 'nsub' = Vnsub)

# Global signal characteristics
TR <- 2
sigmaWhiteNoise <- 100

# Between subject variability: no reasoning for value.
BSubVar <- 2


##
###############
### Monte-Carlo Simulation
###############
##


# For loop over number of parameters
for(j in 1:dim(combPar)[1]){
  # Length of signal (number of scans)
  nscan <- combPar[j,'nscan']
  
  # Number of subjects
  nsub <- combPar[j,'nsub']
  
  # Signal characteristics
  total <- TR*nscan
  on1 <- seq(1,total,40)
  onsets <- list(on1)
  duration <- list(20)
  
  #######################################
  #### DESIGN AND SIGNAL TIME SERIES ####
  #######################################
  
  # Beta_0 and beta_1 (i.e. %BOLD change)
  beta0 <- 100
  BOLDC <- 3
  
  # Generating a design matrix: needed for X
  X <- neuRosim::simprepTemporal(total,1,onsets = onsets,
                                 effectsize = 1, durations = duration,
                                 TR = TR, acc = 0.1, hrf = "double-gamma")
  
  # Generate time series for ONE active voxel: predicted signal, this is the design
  X_value <- neuRosim::simTSfmri(design=X, base=0, SNR=1, noise="none", verbose=FALSE)
  
  # DESIGN PARAMETER NEEDED TO CALCULATE TRUE EFFECT SIZE
  # Extend the design matrix with an intercept
  xIN <- cbind(1,X_value)
  
  # Contrast: not interested in intercept
  CONTRAST <- matrix(c(0,1),nrow=1)
  
  # Calculate (X'X)^(-1) with contrast
  design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)
  
  # Empty vectors for all subjects
  COPE <- GSUB <- VARGSUB <- array(NA,dim=nsub)
  
  # For loop over all subjects
  for(s in 1:nsub){
    # First create subject specific BOLD effect: also depends on number of scans
    signal_BOLDC_subj <- beta0 + 
      (BOLDC + rnorm(1, mean = 0, sd = sqrt(BSubVar))) * X_value
    
    # Generate data for this subject
    Y.subj <- signal_BOLDC_subj + rnorm(n = nscan, mean = 0, sd = sigmaWhiteNoise)
    
    ####************####
    #### ANALYZE DATA: 1e level GLM
    ####************####
    # COPE (beta 1) --> fit GLM
    model.lm <- lm(Y.subj ~ X_value)
    b1 <- coef(model.lm)['X_value']
    COPE[s] <- b1
    
    ####************####
    #### STANDARDIZED ES: 1e level GLM
    ####************####
    # To estimate g, we use the t-value
    Tsub <- summary(model.lm)$coefficients['X_value','t value']
    
    # Estimate g and its variance: exact version
    GSUB[s] <- NeuRRoStat::hedgeG(t = Tsub, N = nscan, type = 'exact')
    VARGSUB[s] <- NeuRRoStat::varHedgeT(g = GSUB[s], N = nscan)
  }
  
  ####************####
  #### STANDARDIZED ES: 2e level GLM
  ####************####
  
  # Weighted average using a fixed effects model
  weights_allSub <- 1/VARGSUB
  WAvg <- sum(weights_allSub * GSUB)/sum(weights_allSub)
  
  # Its variance:
  WAvg_var <- 1/sum(weights_allSub)
  
  # Using random effects MA
  RMA <- rma(as.vector(GSUB), as.vector(VARGSUB), method = 'DL')
  
  
  ####************####
  #### STANDARDIZED GROUP LEVEL ES
  ####************####
  
  # Estimate group level t-value using linear regression (homoscedastic case: OLS)
  TVal <- summary(lm(COPE ~ 1))$coef[,3]
  
  # Some parameters: beta 1 and sample variance
  beta1 <- summary(lm(COPE ~ 1))$coef[,1]
  S2G <- var(COPE)
  # Is identical to: S2G <- (summary(lm(COPE ~ 1))$coef[2]*sqrt(nsub))**2
  
  # Standardized effect size at group stage
  hedge_group <- NeuRRoStat::hedgeG(t = TVal, N = nsub, type = 'exact')
  varHedge_group <- NeuRRoStat::varHedgeT(g = hedge_group, N = nsub)
  

  ####************####
  #### TRUE VALUES
  ####************####
  # 
  # # --> true value for simga^2 at group level: within and between-subject variability
  # trueSigma2G <- sigmaWhiteNoise^2 * design_factor + BSubVar
  # trueG <- BOLDC / sqrt(trueSigma2G) * 
  #   NeuRRoStat::corrH(nsub, type = 'one.sample')
  # trueVarg <- NeuRRoStat::varHedgeT(g = trueG, N = nsub)
  # 
  # # True variance according to Radua
  # trueVargRad <- NeuRRoStat::varHedge(g = trueG, N = nsub)
  # 
  # # True values using the true responses
  # trueG_trueR <- BOLDC / sqrt(BSubVar) * 
  #   NeuRRoStat::corrH(nsub, type = 'one.sample')
  # trueVarg_trueR <- NeuRRoStat::varHedgeT(g = trueG_trueR, N = nsub)
  # 
  # # Save in data frame
  # objects <- data.frame('sim' = rep(K, length(levels(VarHedgeRes$param))),
  #                       'param' = factor(c('beta1', 'S2G', 'hedge', 'varhedge', 
  #                                          'varhedge_radua',
  #                                          'beta1TR', 'S2GTR', 'hedgeTR', 'varhedgeTR')),
  #                       'value' = c(beta1, S2G, hedge, varHedge, varHedgeRad, 
  #                                   beta1_trueR, S2G_trueR, hedge_trueR, varHedge_trueR),
  #                       'nscan' = rep(nscan, length(levels(VarHedgeRes$param))),
  #                       'nsub' = rep(nsub, length(levels(VarHedgeRes$param))),
  #                       'trueValue' = c(BOLDC, trueSigma2G, trueG, trueVarg, trueVargRad,
  #                                       BOLDC, BSubVar, trueG_trueR, trueVarg_trueR))
  # VarHedgeRes <- bind_rows(VarHedgeRes, objects)
}


# Write objects
#saveRDS(VarHedgeRes, file = paste0(wd, '/VarHedgeRes_BSvar_', K,'.rda'))




































####################
#### TITLE:     The effect of time series length on the variance of the standardized effect size.
#### Contents:
####
#### Source Files:
#### First Modified: 20/03/2018
#### Notes:
#################



##
###############
### Notes
###############
##


# My assumption is that the variance of Hedges' g in fMRI data will decrease
# when the length of the time series increases. 
# More particularly, the variance will be underestimated when the time series is not 
# of sufficient lenght. It will only approach the true value when time increases.

# We will check this here.


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
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/BiasVar/Results'
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
VarHedgeRes <- tibble(sim = integer(),
                param = factor(levels = c('beta1', 'S2G', 'hedge', 'varhedge')),
                value = integer(),
                nscan = numeric(),
                nsub = numeric(),
                trueValue = integer(),
                trueValueRad = integer())

checkVar <- tibble(sim = integer(),
                   variance = integer(),
                   nscan = numeric())

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

# Empty vector over all simulations

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
  
  # Now we create the true BOLD signal: depends on number of scans
  signal_BOLDC <- beta0 + BOLDC * X_value
  
  # DESIGN PARAMETER NEEDED TO CALCULATE TRUE EFFECT SIZE
  # Extend the design matrix with an intercept
  xIN <- cbind(1,X_value)
  
  # Contrast: not interested in intercept
  CONTRAST <- matrix(c(0,1),nrow=1)
  
  # Calculate (X'X)^(-1) with contrast
  design_factor <- CONTRAST %*% (solve(t(xIN) %*% xIN )) %*% t(CONTRAST)
  
  # Empty vectors for all subjects
  COPE <- VARCOPE <- array(NA,dim=nsub)
  
  # For loop over all subjects
  for(s in 1:nsub){
    # Generate data for this subject
    Y.subj <- signal_BOLDC + rnorm(n = nscan, mean = 0, sd = sigmaWhiteNoise)

    ####************####
    #### ANALYZE DATA: 1e level GLM
    ####************####
    # COPE (beta 1) --> fit GLM
    model.lm <- lm(Y.subj ~ X_value)
    b1 <- coef(model.lm)['X_value']
    COPE[s] <- b1
  }
  
  # Estimate group level t-value using linear regression (homoscedastic case: OLS)
  TVal <- summary(lm(COPE ~ 1))$coef[,3]
  
  # Some parameters: beta 1 and sample variance
  beta1 <- summary(lm(COPE ~ 1))$coef[,1]
  S2G <- var(COPE)
      # Is identical to: S2G <- (summary(lm(COPE ~ 1))$coef[2]*sqrt(nsub))**2

  # Standardized effect size
  hedge <- NeuRRoStat::hedgeG(t = TVal, N = nsub, type = 'exact')
  varHedge <- NeuRRoStat::varHedgeT(g = hedge, N = nsub)
  
  # True values
  # --> true value for simga^2 at group level: no between-study variability
  trueSigma2G <- sigmaWhiteNoise^2 * design_factor
  trueG <- BOLDC / (sigmaWhiteNoise * sqrt(design_factor)) * NeuRRoStat::corrH(nsub, type = 'one.sample')
  trueVarg <- NeuRRoStat::varHedgeT(g = trueG, N = nsub)
  
  # True variance according to Radua
  trueVargRad <- NeuRRoStat::varHedge(g = trueG, N = nsub)
  
  # Save in data frame
  objects <- data.frame('sim' = rep(K, length(levels(VarHedgeRes$param))),
               'param' = factor(c('beta1', 'S2G', 'hedge', 'varhedge')),
               'value' = c(beta1, S2G, hedge, varHedge),
               'nscan' = rep(nscan, length(levels(VarHedgeRes$param))),
               'nsub' = rep(nsub, length(levels(VarHedgeRes$param))),
               'trueValue' = c(BOLDC, trueSigma2G, trueG, trueVarg),
               'trueValueRad' = c(NA, NA, NA, trueVargRad))
  VarHedgeRes <- bind_rows(VarHedgeRes, objects)
}


# Write objects
saveRDS(VarHedgeRes, file = paste0(wd, '/VarHedgeRes_', K,'.rda'))

# VarHedgeRes %>% select(-nsub) %>%
#   group_by(param, nscan, trueValue) %>%
#   summarise(avg = mean(value)) %>%
#   ungroup() %>% filter(param == 'hedge') %>%
#   ggplot(., aes(x = nscan, y = avg)) +
#   geom_line() + 
#   geom_line(aes(x = nscan, y = trueValue), colour = 'green') +
#   scale_y_continuous('Var(g)')
# 
# 
# checkVar %>% group_by(nscan) %>%
#   summarise(avgVar = mean(variance)) %>%
#   ggplot(., aes(x = nscan, y = avgVar)) +
#   geom_line()





































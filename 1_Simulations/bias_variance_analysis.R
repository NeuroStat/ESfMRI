####################
#### TITLE:     Analysis of the effect of time series length and sample size
####            on the variance of the standardized effect size.
#### Contents:
####
#### Source Files:
#### First Modified: 22/03/2018
#### Notes:
#################



##
###############
### Notes
###############
##


# Variance of g increases as g increases. Which is the case when the time
# series is larger.
# But this effect might change when the sample size increases.
# We will check this here.


##
###############
### Preparation
###############
##


# Provide data location
DataLoc <- '~/Desktop/Results/'


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

# Empty data frame with results
VarHedgeRes <- tibble(sim = integer(),
                      param = factor(levels = c('beta1', 'S2G', 'hedge', 'varhedge')),
                      value = integer(),
                      nscan = numeric(),
                      nsub = numeric(),
                      trueValue = integer(),
                      trueValueRad = integer())

# Number of simulations
nsim <- 1000


##
###############
### Read in data
###############
##

# For loop over the simulations
for(i in 1:nsim){
  # Read in data object and bind to data frame
  VarHedgeRes <- bind_rows(VarHedgeRes,
    readRDS(paste0(DataLoc, 'VarHedgeRes_', i, '.rda')))
}


##
###############
### Plot results
###############
##


# Average over simulations
avgRes <- VarHedgeRes %>% group_by(param, nscan,
                            nsub, trueValue, trueValueRad) %>%
    summarise(AvgValue = mean(value)) %>% ungroup()

# Plot beta
avgRes %>% filter(param == 'beta1') %>%
  ggplot(., aes(x = nsub, y = AvgValue, group = nscan)) +
  geom_line(aes(colour = nscan))


# Boxplots of beta estimates
VarHedgeRes %>% filter(param == 'beta1') %>%
  filter(nscan %in% c(100,130,260,500)) %>%
  ggplot(., aes(x = factor(nsub), y = value)) +
  geom_boxplot(outlier.size = .5) +
  geom_hline(aes(yintercept = 
                   VarHedgeRes %>% filter(param == 'beta1') %>% 
                   select(trueValue) %>%
                   unique() %>% as.numeric()),
             linetype = 'dotted') +
  scale_x_discrete('Number of subjects') +
  scale_y_continuous(expression(hat(beta))) +
  facet_grid(nscan ~ .) +
  labs(subtitle = "Dotted line represents true value") +
  theme_bw()


# Sample variance
avgRes %>% filter(param == 'S2G') %>%
  filter(nsub %in% c(20,50,70,100)) %>%
  ggplot(., aes(x = nscan, y = AvgValue)) +
  geom_line() +
  facet_grid(nsub ~ .)
  geom_line(aes(colour = nsub))
  
# Boxplots of variance
VarHedgeRes %>% filter(param == 'S2G') %>%
  filter(nsub %in% c(20,50,70,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_boxplot(outlier.size = .5) +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param == 'S2G') %>% 
                 select(nscan, nsub, trueValue) %>%
                 filter(nsub %in% c(20,50,70,100)) %>%
                 filter(nscan %in% seq(100,500,by = 40)) %>%
                 distinct(),
               aes(x = nscan, y = trueValue), colour = 'orange') +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(S[G]^2)) +
  facet_grid(nsub ~ .) +
  labs(subtitle = "Dotted line represents true value") +
  theme_bw()

VarHedgeRes %>% filter(param == 'S2G') %>% 
  select(nscan, nsub, trueValue) %>%
  filter(nsub %in% c(20,50,70,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  distinct()
 

unique(VarHedgeRes$nscan)
  
seq(100,500,by = 40)








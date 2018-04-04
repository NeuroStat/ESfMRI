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
DataLoc_list <- list('BSvar0' = "/Volumes/Elements/ReviewESfMRI/Simulations/BSVar0",
                "BSVar2" = "/Volumes/Elements/ReviewESfMRI/Simulations/BSVar2")
DataLoc_list <- list('BSvar0' = "/Volumes/Elements/ReviewESfMRI/Simulations/BSVar0",
                     "BSVar2" = "~/Desktop/ResultsB",
                     "BSVar3" = "~/Desktop/Results")
dataWD <- "BSVar3"
DataLoc <- DataLoc_list[[dataWD]]


# Load in libraries
library(lattice)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(reshape2)
library(cowplot)
library(ggridges)
library(RColorBrewer)
library(mvmeta)
library(metafor)
library(devtools)
library(neuRosim)
library(NeuRRoStat)

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


##
###############
### Functions
###############
##

# Function to plot Ridge plots with true values
    # data are the simulation results
    # parameter: which parameter do you want to plot?
    # N: vector of sample sizes
    # xAxis: label for x-axis
    # offSet: vector of two offset numbers to add to x-axis if needed
RidgeTrueV <- function(data, parameter, N, xAxis, offSet = 0){
  # Check length of N
  if(length(N) != 2) stop('Function created for two sample sizes')
  
  # We need a character vector of sample sizes as well (for plot)
  Nchar <- as.character(N)
  
  # First create data frame with the true values
  TrueData <- VarHedgeRes %>% filter(param == parameter) %>% 
    select(nscan, nsub, trueValue) %>%
    filter(nsub %in% N) %>%
    filter(nscan %in% seq(100, 500, by = 40)) %>%
    distinct()
  TrueData$nscan <- factor(TrueData$nscan)
  
  # Here I select the minimum and 99th quantile of the data to cut-off the axes.
  # Otherwise the tails skew the visualization.
  xLims <- VarHedgeRes %>% filter(param == parameter) %>%
    filter(nsub %in% N) %>%
    filter(nscan %in% seq(100,500,by = 40)) %>% 
    select(value) %>% summarise(min_value = min(value),
                                perc99_value = quantile(value, p = 0.99))
  
  # Adding the offset values to the xLims
  xLims <- xLims + offSet
  
  # Colour scheme according to subjects
  subCol <- data.frame(N = seq(20, 100, by = 10),
                       colour = 
    c('#41b6c4','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
      '#045a8d','#3690c0','#02818a','#023858'),
    stringsAsFactors = FALSE)
  col1 <- subCol[which(subCol[,'N'] == N[1]), 'colour']
  col2 <- subCol[which(subCol[,'N'] == N[2]), 'colour']

  # Now we create the Ridge plot
  ToPlot <- VarHedgeRes %>% filter(param == parameter) %>%
    filter(nsub %in% N) %>%
    filter(nscan %in% seq(100,500,by = 40)) %>%
    ggplot(., aes(x = value, y = factor(nscan), fill = factor(nsub))) +
    geom_density_ridges(scale=0.9, color='white') +
    scale_fill_cyclical(breaks = Nchar,
                        values = c(col1, col2),
                        name = "N", guide = "legend") +
    geom_segment(data = TrueData,
                 aes(x = trueValue, xend = trueValue, 
                     y = as.numeric(nscan), yend = as.numeric(nscan) + 0.9),
                 color = "red") +
    scale_x_continuous(xAxis, limits = c(as.numeric(xLims['min_value']), 
                                            as.numeric(xLims['perc99_value']))) +
    scale_y_discrete("Number of scans", expand = c(0.01, 0)) 
  
  return(ToPlot)
}

# Ridge plots using N as facets
RidgeTrueV_facet <- function(data, parameter, N, xAxis){
  # Check length of N
  if(length(N) != 2) stop('Function created for two sample sizes')
  
  # We need a character vector of sample sizes as well (for plot)
  Nchar <- as.character(N)
  
  # Data frame with true values 
  TrueData <- VarHedgeRes %>% filter(param == parameter) %>% 
    select(nscan, nsub, trueValue) %>%
    filter(nsub %in% N) %>%
    filter(nscan %in% seq(100, 500, by = 40)) %>%
    distinct()
  TrueData$nscan <- factor(TrueData$nscan)
  
  # Colour scheme according to subjects
  subCol <- data.frame(N = seq(20, 100, by = 10),
                       colour = 
                         c('#41b6c4','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                           '#045a8d','#3690c0','#02818a','#023858'),
                       stringsAsFactors = FALSE)
  col1 <- subCol[which(subCol[,'N'] == N[1]), 'colour']
  col2 <- subCol[which(subCol[,'N'] == N[2]), 'colour']
  
  # Adapt data frame
  DataToPlot <- VarHedgeRes %>% filter(param == parameter) %>%
    filter(nsub %in% N) %>%
    filter(nscan %in% seq(100, 500,by = 40)) %>%
    mutate(nsubLabel = paste('N = ', nsub, sep = '')) 
  DataToPlot$nsubLabel <- factor(DataToPlot$nsubLabel, 
                                 levels = paste0('N = ', N))
  
  # Now create the facet plot
  ToPlot <- DataToPlot %>%
    ggplot(., aes(x = value, y = factor(nscan), fill = factor(nsub))) +
    geom_density_ridges(scale=0.9, color='white') +
    geom_segment(data = TrueData,
                 aes(x = trueValue, xend = trueValue, 
                     y = as.numeric(nscan), yend = as.numeric(nscan) + 0.9),
                 color = "red") +
    facet_grid(. ~ nsubLabel) +
    scale_fill_manual(values = c(col1, col2)) +
    scale_x_continuous(xAxis) +
    scale_y_discrete("Number of scans", expand = c(0.01, 0)) +
    theme(legend.position = "none")
  
  return(ToPlot)
}


##
###############
### Read in data
###############
##

# For loop over the simulations
for(i in 1:nsim){
  # Read in data object and bind to data frame
  VarHedgeRes <- bind_rows(VarHedgeRes,
    readRDS(paste0(DataLoc, '/VarHedgeRes_BSvar_', i, '.rda')))
}

# Average over simulations
avgRes <- VarHedgeRes %>% group_by(param, nscan,
                           nsub, trueValue) %>%
  summarise(AvgValue = mean(value)) %>% ungroup()


##
###############
### Plot results
###############
##


# Plot beta: two-level approach
avgRes %>% filter(param == 'beta1') %>%
  ggplot(., aes(x = nsub, y = AvgValue, group = nscan)) +
  geom_line(aes(colour = nscan))


# Boxplots of beta estimates
VarHedgeRes %>% filter(param == 'beta1') %>%
  filter(nscan %in% c(100,130,260,500)) %>%
  ggplot(., aes(x = factor(nsub), y = value)) +
  geom_boxplot(outlier.size = .5,
               position = 'dodge2') +
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

# Boxplots of beta estimates according to model
VarHedgeRes %>% filter(param %in% c('beta1', 'beta1TR')) %>%
  filter(nscan %in% c(100,260,500)) %>%
  ggplot(., aes(x = factor(nsub), y = value, fill = param)) +
  geom_boxplot(aes(fill = param), outlier.size = .5,
               position = 'dodge') +
  geom_hline(aes(yintercept = 
                   VarHedgeRes %>% filter(param == 'beta1') %>% 
                   select(trueValue) %>%
                   unique() %>% as.numeric()),
             linetype = 'dotted') +
  scale_x_discrete('Number of subjects') +
  scale_y_continuous(expression(hat(beta))) +
  scale_fill_manual('Fitted model', 
                    values = c('#f4a582', '#2166ac'),
                    labels = c('Two-stage', 'True responses')) +
  facet_grid(nscan ~ .) +
  labs(subtitle = "Dotted line represents true value") +
  theme_bw() +
  theme(legend.position = 'bottom')

# Violin of beta estimates according to model
VarHedgeRes %>% filter(param %in% c('beta1', 'beta1TR')) %>%
  filter(nscan %in% c(100,260,500)) %>%
  mutate(paramLabel = recode(param,
                          beta1 = 'Two-stage model',
                          beta1TR = 'True responses'),
         nscanLabel = recode(nscan,
                             '100' = 'T = 100',
                             '260' = 'T = 260',
                             '500' = 'T = 500')) %>%
  ggplot(., aes(x = factor(nsub), y = value)) +
  geom_violin() +
  geom_hline(aes(yintercept = 
                   VarHedgeRes %>% filter(param == 'beta1') %>% 
                   select(trueValue) %>%
                   unique() %>% as.numeric()),
             linetype = 'dotted') +
  scale_x_discrete('Number of subjects') +
  scale_y_continuous(expression(hat(beta))) +
  facet_grid(nscanLabel ~ paramLabel) +
  labs(subtitle = "Dotted line represents true value") +
  theme_bw() +
  theme(legend.position = 'bottom')


# Sample variance
avgRes %>% filter(param == 'S2G') %>%
  filter(nsub %in% c(20,50,70,100)) %>%
  ggplot(., aes(x = nscan, y = AvgValue)) +
  geom_line() +
  facet_grid(nsub ~ .)

# Boxplots of variance
VarHedgeRes %>% filter(param == 'S2G') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_boxplot(outlier.size = .5) +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param == 'S2G') %>% 
                 select(nscan, nsub, trueValue) %>%
                 filter(nsub %in% c(20,100)) %>%
                 filter(nscan %in% seq(100,500,by = 40)) %>%
                 distinct(),
               aes(x = factor(nscan), y = trueValue), colour = 'orange', size = 0.3) +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(S[G]^2)) +
  facet_grid(nsub ~ .) +
  labs(subtitle = "Solid orange line represents true value") +
  theme_bw()

# Ridge plots for N = 20 and 100, using facets
xAxis <- expression(S[G]^2)
RidgeTrueV_facet(data = VarHedgeRes, parameter = 'S2G', N = c(20, 100),
                 xAxis = xAxis)

# Estimate of g
avgRes %>% filter(param == 'hedge') %>%
  filter(nsub %in% c(20,50,70,100)) %>%
  ggplot(., aes(x = nscan, y = AvgValue)) +
  geom_line() +
  facet_grid(nsub ~ .)


# Check whether g and var(g) are unbiased
VarHedgeRes %>% filter(param %in% c('hedge', 'varhedge') ) %>%
  filter(nsub %in% c(20,50,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  group_by(nsub, nscan, param, trueValue) %>%
  summarise(AvgEmp = mean(value)) %>% 
  tidyr::gather(., key = "type", value = "value", 4:5) %>%
  ungroup() %>% 
  mutate(LabelNsub = paste('N = ', nsub, sep = ''),
         LabelParam = recode_factor(param,
           'hedge' = "Hedges' g",
           'varhedge' = 'Var(g)')) %>%
  group_by(type) %>%
  ggplot(., aes(x = nscan, y = value)) +
  geom_col(aes(fill = type), position = 'dodge') +
  facet_wrap(LabelNsub ~ LabelParam, scales = 'free', ncol = 2) +
  scale_fill_manual("", values = c('#b3e2cd', '#fdcdac'),
                    label = c('Average over Monte-Carlo simulations',
                              'Expected value')) +
  scale_x_continuous("Number of scans") + 
  scale_y_continuous("") +
  theme_bw() +
  theme(legend.position = 'top')

# Create a density of Hedges g with the expected density (t-distribution) 
    # for N = 100 and nscan = 500 versus the observed distribution
TrueN100 <- VarHedgeRes %>% filter(param == 'hedge') %>%
  filter(nsub == 100 ) %>%
  filter(nscan == 500) %>%
  select(trueValue) %>%
  distinct()
expDensity <- 1/sqrt(100) * rt(n = 1000, df = 99, ncp = as.numeric(sqrt(100)*TrueN100))
VarHedgeRes %>% filter(param == 'hedge') %>%
  filter(nsub == 100 ) %>%
  filter(nscan == 500) %>%
  mutate(expDensity = expDensity) %>%
  select(sim, value, expDensity) %>%
  rename(obsDensity = value) %>%
  tidyr::gather(., key = 'type', value = 'value', 2:3) %>%
  ggplot(., aes(x = value, colour = type)) +
  geom_density(size = 1.2) +
  scale_color_brewer("Density", labels = c('Expected', 'Observed'),
                     type = 'qual', palette = 2)

# Boxplot of Hedges estimates
VarHedgeRes %>% filter(param == 'hedge') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  mutate(nsubLabel = recode_factor(nsub,
                       '20' = 'N = 20',
                       '100' = 'N = 100')) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_boxplot(outlier.size = .5) +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param == 'hedge') %>% 
                 select(nscan, nsub, trueValue) %>%
                 filter(nsub %in% c(20,100)) %>%
                 filter(nscan %in% seq(100,500,by = 40)) %>%
                 distinct(),
               aes(x = factor(nscan), y = trueValue), colour = 'orange', size = 0.3) +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(g[e])) +
  facet_grid(nsubLabel ~ .) +
  labs(subtitle = "Solid orange line represents true value") +
  theme_bw()

# Or using ridge plots
xAxis <- expression(g[e])
RidgeTrueV_facet(data = VarHedgeRes, parameter = 'hedge', N = c(20, 100),
                 xAxis = xAxis)

# Violin plot of variance of g according to subjects
VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_violin(scale = 'area') + 
  facet_grid(nsub ~ ., scales = 'free_y')

# Boxplot of variance of g
VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_boxplot(outlier.size = .5) +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param == 'varhedge') %>% 
                 select(nscan, nsub, trueValue) %>%
                 filter(nsub %in% c(20,100)) %>%
                 filter(nscan %in% seq(100, 500, by = 40)) %>%
                 distinct(),
               aes(x = factor(nscan), y = trueValue), colour = 'orange', size = 0.3) +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(Var(g[e]))) +
  facet_grid(nsub ~ ., scales = 'free_y') +
  labs(subtitle = "Solid orange line represents true value") +
  theme_bw()

# Same but using the formula of Radua
VarHedgeRes %>% filter(param == 'varhedge_radua') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_boxplot(outlier.size = .5) +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param == 'varhedge_radua') %>% 
                 select(nscan, nsub, trueValue) %>%
                 filter(nsub %in% c(20,100)) %>%
                 filter(nscan %in% seq(100, 500, by = 40)) %>%
                 distinct(),
               aes(x = factor(nscan), y = trueValue), colour = 'orange', size = 0.3) +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(Var(g[e]))) +
  facet_grid(nsub ~ ., scales = 'free_y') +
  labs(subtitle = "Solid orange line represents true value") +
  theme_bw()

# Now using the true responses
VarHedgeRes %>% filter(param == 'varhedgeTR') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_boxplot(outlier.size = .5) +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param == 'varhedgeTR') %>% 
                 select(nscan, nsub, trueValue) %>%
                 filter(nsub %in% c(20,100)) %>%
                 filter(nscan %in% seq(100, 500, by = 40)) %>%
                 distinct(),
               aes(x = factor(nscan), y = trueValue), colour = 'orange', size = 0.3) +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(Var(g[e]))) +
  facet_grid(nsub ~ ., scales = 'free_y') +
  labs(subtitle = "Solid orange line represents true value") +
  theme_bw()

# Fitted model and true responses
VarHedgeRes %>% filter(param %in% c('varhedge', 'varhedgeTR')) %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value, fill = param)) +
  geom_boxplot(outlier.size = .5) +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param %in% c('varhedge', 'varhedgeTR')) %>% 
                 select(nscan, nsub, trueValue, param) %>%
                 filter(nsub %in% c(20,100)) %>%
                 filter(nscan %in% seq(100, 500, by = 40)) %>%
                 distinct(),
               aes(x = factor(nscan), y = trueValue, colour = param), size = 0.3) +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(Var(g[e]))) +
  facet_grid(nsub ~ ., scales = 'free_y') +
  labs(subtitle = "Solid coloured line represents true value") +
  theme_bw()

# Violin plot of variance of g according to model
VarHedgeRes %>% filter(param %in% c('varhedge', 'varhedgeTR')) %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  mutate(paramLabel = recode(param,
                             varhedge = 'Two-stage model',
                             varhedgeTR = 'True responses'),
         nsubLabel = recode(nsub,
                             '20' = 'N = 20',
                             '100' = 'N = 100')) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_violin() +
  scale_x_discrete('Number of subjects') +
  scale_y_continuous(expression(hat(beta))) +
  facet_wrap(nsubLabel ~ paramLabel, ncol = 2, scales = 'free') +
  labs(subtitle = "Dotted line represents true value") +
  theme_bw() +
  theme(legend.position = 'bottom')


# Violin plot of variance of g
VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_violin() +
  geom_boxplot(data = 
                 VarHedgeRes %>% filter(param == 'varhedge') %>% 
                 select(nscan, nsub, trueValue) %>%
                 filter(nsub %in% c(20,100)) %>%
                 filter(nscan %in% seq(100, 500, by = 40)) %>%
                 distinct(),
               aes(x = factor(nscan), y = trueValue), colour = 'orange', size = 0.3) +
  scale_x_discrete('Number of scans') +
  scale_y_continuous(expression(Var(g[e]))) +
  facet_grid(nsub ~ ., scales = 'free_y') +
  labs(subtitle = "Solid orange line represents true value") +
  theme_bw()


# Ridge plot on variance with N = 90 and 100
# Have some offSet for left side of axis. 
offSet <- c(-0.0003, 0)
secL_varG <- RidgeTrueV(data = VarHedgeRes, parameter = "varhedge", N = c(90,100), 
           xAxis = 'Var(g)', offSet = offSet)

# Using Radua his formula
RidgeTrueV(data = VarHedgeRes, parameter = "varhedge_radua", N = c(90,100), 
           xAxis = 'Var(g)', offSet = offSet)

# On true responses with N = 50 and 100 (90 is too close to 100 for visualization)
offSet <- c(0, 0)
TR_varG <- RidgeTrueV(data = VarHedgeRes, parameter = "varhedgeTR", N = c(50,100), 
           xAxis = 'Var(g)', offSet = offSet)

# Have both ridge plots in one figure
plot_grid(secL_varG, TR_varG, labels = c("A", "B"), nrow = 2, align = "v")







# Older code

data.frame(x = 1/sqrt(100) * rt(n = 1000, df = 99, ncp = sqrt(100)*0.292)) %>%
  ggplot(., aes(x = x)) +
  geom_density()

VarHedgeRes %>% filter(param == 'hedge') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan == 500) %>%
  group_by(nsub) %>%
  summarise(var(value))
# Variantie
VarHedgeRes %>% filter(param == 'hedge') %>%
  filter(nsub == 100 ) %>%
  filter(nscan == 500) %>%
  summarise(var(value))

VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub == 100 ) %>%
  filter(nscan == 500) %>%
  summarise(mean(value))

VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub == 100 ) %>%
  filter(nscan == 500) %>%
  select(trueValue) %>%
  distinct()


# Gemiddelde
VarHedgeRes %>% filter(param == 'hedge') %>%
  filter(nsub == 100 ) %>%
  filter(nscan == 500) %>%
  summarise(mean(value))

VarHedgeRes %>% filter(param == 'hedge') %>% 
  select(nscan, nsub, trueValue) %>%
  filter(nsub == 100) %>%
  filter(nscan == 500 ) %>%
  distinct()

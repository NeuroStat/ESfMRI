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
                     "BSVar3" = "~/Desktop/")
dataWD <- "BSVar2"
DataLoc <- DataLoc_list[[dataWD]]


# Load in libraries
library(lattice)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tibble)
library(reshape2)
library(ggridges)
library(RColorBrewer)
library(mvmeta)
library(metafor)
library(devtools)
library(neuRosim)
library(NeuRRoStat)

# Empty data frame with results
VarHedgeRes <- tibble(sim = integer(),
                      param = factor(levels = c('beta1', 'S2G', 'hedge', 'varhedge', 
                              'varhedge_radua',
                              'beta1TR', 'S2GTR', 'hedgeTR', 'varhedgeTR')),
                      value = integer(),
                      nscan = numeric(),
                      nsub = numeric(),
                      trueValue = integer())

# Number of simulations
nsim <- 500


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
VarHedgeRes %>% filter(param %in% c('beta1', 'beta1TR')) %>%
  filter(nscan %in% c(100,130,260,500)) %>%
  ggplot(., aes(x = factor(nsub), y = value, fill = param)) +
  geom_boxplot(aes(fill = param), outlier.size = .5,
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


# Estimate of g
avgRes %>% filter(param == 'hedge') %>%
  filter(nsub %in% c(20,50,70,100)) %>%
  ggplot(., aes(x = nscan, y = AvgValue)) +
  geom_line() +
  facet_grid(nsub ~ .)

# Boxplot of Hedges estimates
VarHedgeRes %>% filter(param == 'hedge') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
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
  facet_grid(nsub ~ .) +
  labs(subtitle = "Solid orange line represents true value") +
  theme_bw()


# Variance of g
avgRes %>% filter(param == 'varhedge') %>%
  filter(nsub %in% c(20,50,70,100)) %>%
  ggplot(., aes(x = nscan, y = AvgValue)) +
  geom_line() +
  facet_grid(nsub ~ .)

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


# Same but using the true responses
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



# Try other figures
VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = factor(nscan), y = value)) +
  geom_violin(scale = 'area') + 
  facet_grid(nsub ~ ., scales = 'free_y')

VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub == 100) %>%
  #filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = value, y = factor(nscan))) +
  geom_density_ridges(scale = 0.9) + theme_ridges() +
  scale_x_continuous(expression(Var(g[e])), expand = c(0.01, 0)) +
  scale_y_discrete('Number of scans', expand = c(0.01, 0)) +
  geom_segment(data = VarHedgeRes %>% filter(param == 'varhedge') %>% 
                 select(nscan, nsub, trueValue) %>%
                 #filter(nsub %in% c(20,100)) %>%
                 filter(nsub == 100) %>%
                 filter(nscan %in% seq(100, 500, by = 40)) %>%
                 distinct(),
                 aes(x = trueValue, xend = trueValue, 
                   y = as.numeric(nscan), yend = as.numeric(nscan) + 0.9)) + 
  #geom_vline(aes(xintercept = trueValue)) + 
  ggtitle('Density of estimated variance of g')
  


iris_lines <- data.frame(Species = c("setosa", "versicolor", "virginica"),
                         x0 = c(5, 4, 7))

ggplot(iris, aes(x=Sepal.Length, y=Species, fill=..x..)) +
  ?geom_density_ridges_gradient(jittered_points = FALSE, quantile_lines = 
                                 FALSE, quantiles = 2, scale=0.9, color='white') +
  geom_segment(data = iris_lines, aes(x = x0, xend = x0, y = as.numeric(Species),
                                      yend = as.numeric(Species) + .9),
               color = "red") +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme_ridges(grid = FALSE, center = TRUE)


TrueData <- VarHedgeRes %>% filter(param == 'varhedge') %>% 
  select(nscan, nsub, trueValue) %>%
  #filter(nsub %in% c(20,100)) %>%
  filter(nsub == 100) %>%
  filter(nscan %in% seq(100, 500, by = 40)) %>%
  distinct()
TrueData$nscan <- factor(TrueData$nscan)


VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub == 20) %>%
  #filter(nsub %in% c(20,100)) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
ggplot(., aes(x = value, y = factor(nscan))) +
  geom_density_ridges(scale=0.9, color='white') +
  geom_segment(data = TrueData,
               aes(x = trueValue, xend = trueValue, 
                   y = as.numeric(nscan), yend = as.numeric(nscan) + 0.9),
               color = "red") +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme_ridges(grid = FALSE, center = TRUE)


VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub == 100) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = value, y = factor(nscan))) +
  geom_density_ridges(scale=0.9, color='white') +
  geom_segment(data = TrueData,
               aes(x = trueValue, xend = trueValue,
                   y = as.numeric(nscan), yend = as.numeric(nscan) + 0.9),
               color = "red") +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme_ridges(grid = FALSE, center = TRUE)

# On true responses
TrueData <- VarHedgeRes %>% filter(param == 'varhedgeTR') %>% 
  select(nscan, nsub, trueValue) %>%
  #filter(nsub %in% c(20,100)) %>%
  filter(nsub == 100) %>%
  filter(nscan %in% seq(100, 500, by = 40)) %>%
  distinct()
TrueData$nscan <- factor(TrueData$nscan)

VarHedgeRes %>% filter(param == 'varhedgeTR') %>%
  filter(nsub == 100) %>%
  filter(nscan %in% seq(100,500,by = 40)) %>%
  ggplot(., aes(x = value, y = factor(nscan))) +
  geom_density_ridges(scale=0.9, color='white') +
  geom_segment(data = TrueData,
               aes(x = trueValue, xend = trueValue,
                   y = as.numeric(nscan), yend = as.numeric(nscan) + 0.9),
               color = "red") +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme_ridges(grid = FALSE, center = TRUE)

# Repeat for hedges g



VarHedgeRes %>% filter(param == 'varhedge') %>%
  filter(nsub == 20) %>%
  filter(nscan  == 100) %>%
  ggplot(., aes(x = value)) + 
  geom_density() +
  geom_vline(data = 
               VarHedgeRes %>% filter(param == 'varhedge') %>%
               filter(nsub == 20) %>%
               filter(nscan  == 100) %>%
               select(trueValue) %>%
               unique(),
             aes(xintercept = trueValue))
  



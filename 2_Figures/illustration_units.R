####################
#### TITLE: Combine fMRI time series with different units.
#### Contents: 	
#### 
#### Source Files: 
#### First Modified: 21/11/2018
#### Notes: 
#################

# Illustrations about combining fMRI time series with
# differences in units.

##
###############
### Preparation
###############
##

# Load in libraries
library(tidyverse)
library(cowplot)
library(neuRosim)
library(NeuRRoStat)

# Set seed
set.seed(111290)

# Function to generate data (time series): scaled according to percent BOLD change
genDat <- function(base, BOLDperc, X, nscan, sigmaWN){
  series <- base + (base * (BOLDperc/100) * (X + rnorm(n = nscan, sd = sigmaWN)))  
  return(series)
}




##
###############
### Parameters
###############
##

# Number of simulations
nsim <- 1000

# Different base values for actual data generation
base1 <- 10000
base2 <- 100

# Number of subjects
N <- 50

# Signal characteristics
TR <- 2
nscan <- 200
total <- TR*nscan
on1 <- seq(1,total,40)
onsets <- list(on1)
duration <- list(20)
BOLDes <- 1

# True BOLD response in %
BOLDperc <- 3

# Generating a design matrix --> preparation
Xprep <- neuRosim::simprepTemporal(totaltime = total,
       regions = 1,
       onsets = onsets,
       effectsize = BOLDes, 
       durations = duration,
       TR = TR, acc = 0.1, 
       hrf="double-gamma")

# Generate the actual block design time series 
X <- neuRosim::simTSfmri(design=Xprep,
                         base = 0,
                         SNR = 1, 
                         noise="none", verbose=FALSE)

plot(X, type = 'l')

# Noise characteristics (within and between subject resp.)
sigmaWN <- 2
bsSD <- 1


##
###############
### Illustration of single subjects
###############
##

# First for illustration: 3 time series.
# Both signals have the same BOLD % signal change.
# So we scale according to baseline.
# Note that we scale the sigma as well!
illu1 <- genDat(base = base1, BOLDperc = BOLDperc, X = X, 
                nscan = nscan, sigmaWN = sigmaWN)
illu2 <- genDat(base = base2, BOLDperc = BOLDperc, X = X, 
                nscan = nscan, sigmaWN = sigmaWN)

# Let us plot the observed time series with the block events
par(mfrow = c(1,3))
plot(illu1, type = 'l', xlab = 'scans', ylab = '', main = 'Normalization 1')
lines((base1 + (base1 * (BOLDperc/100) * X)), col = 'green', 
      lty = 'dashed', lwd = 2)
plot(illu2, type = 'l', xlab = 'scans', ylab = '', main = 'Normalization 2')
lines((base2 + (base2 * (BOLDperc/100) * X)), col = 'green', 
      lty = 'dashed', lwd = 2)
plot(illu1, ylim = c(min(illu1, illu2), max(illu1, illu2)),
     type = 'l', xlab = 'scans', ylab = '', main = 'Both normalizations')
lines(illu2)
par(mfrow = c(1,1))

# Let us create a data frame to use with ggplot
visTimS <- data.frame('timeSeries' = c(illu1,
                          (base1 + (base1 * (BOLDperc/100) * X)),
                          illu2,
                          (base2 + (base2 * (BOLDperc/100) * X))),
           'type' = rep(c('obsResp', 'trueResp'), each = nscan),
           'norm' = rep(c(1,2), each = nscan * 2),
           'scan' = rep(c(1:nscan),4))


# Create 3 plots and combine them with gggrid (too cumbersome with facet wrap)
norm1 <- visTimS %>%
  filter(norm == 1) %>%
  mutate(typeF = factor(type)) %>%
  ggplot(., aes(x = scan, y = timeSeries, group = typeF)) + 
  geom_line(aes(colour = typeF, linetype = typeF, size = typeF)) + 
  scale_color_manual('', values = c('black', '#3690c0'),
                     labels = c('observed response',
                                'true response')) +
  scale_size_manual('', values = c(0.3,1.2),
                    labels = c('observed response',
                               'true response')) +
  scale_linetype_manual('', values = c('solid', 'dashed'),
                        labels = c('observed response',
                                   'true response')) +
  ggtitle('normalization 1') +
  scale_y_continuous('BOLD response') + 
  scale_x_continuous('scans') +
  guides(size = FALSE) +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
norm1

norm2 <- visTimS %>%
  filter(norm == 2) %>%
  mutate(typeF = factor(type)) %>%
  ggplot(., aes(x = scan, y = timeSeries, group = typeF)) + 
  geom_line(aes(colour = typeF, linetype = typeF, size = typeF)) + 
  scale_color_manual('', values = c('black', '#3690c0'),
                     labels = c('observed response',
                                'true response')) +
  scale_size_manual('', values = c(0.3,1.2),
                    labels = c('observed response',
                               'true response')) +
  scale_linetype_manual('', values = c('solid', 'dashed'),
                        labels = c('observed response',
                                   'true response')) +
  ggtitle('normalization 2') +
  scale_y_continuous('BOLD response') + 
  scale_x_continuous('scans') +
  guides(size = FALSE) +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
norm2

# Both on transformed scale
norm3 <- visTimS %>%
  filter(type == 'obsResp') %>%
  ggplot(., aes(x = scan, y = timeSeries, group = norm)) + 
  geom_line() + 
  ggtitle('both normalizations') +
  scale_y_continuous('log(BOLD) response') + 
  scale_x_continuous('scans') +
  coord_trans(y = "log10") +
  guides(size = FALSE) +
  labs(caption = 'Each signal has a 3% BOLD signal change') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
norm3

# Now group the 3 panels with one shared legend
legend_b <- get_legend(norm1 + theme(legend.position="bottom"))
grid3P <- plot_grid(norm1 + theme(legend.position="none"), 
                    norm2 + theme(legend.position="none"), 
                    norm3 + theme(legend.position="none"),
                    nrow = 1, align = "hv",
                    axis = 'tblr')
plot_grid(grid3P, legend_b, ncol = 1, rel_heights = c(1, .1))

# Save plot
ggsave(filename = paste(getwd(), '/illu_norm.png', sep = ''), 
       plot = last_plot(), 
       dpi = 320, scale = 1.5, 
       width = 16, height = 8, units = 'cm')

##
###############
### Generate data: separate normalizations
###############
##

# Generate data for all subjects and then combine them using OLS
# Small amount of between subject variability!

# Empty data frame over all simulations
dataGroup <-
  data.frame() %>% as_tibble()

# For loop over the amount of simulations
for(k in 1:nsim){
  # Empty data frame for subjects
  dataSubj <- data.frame() %>% as_tibble()
  
  # For loop over the two normalizations
  for(j in 1:2){
    if(j == 1) baseJ <- base1
    if(j == 2) baseJ <- base2
    # For loop over all subjects
    for(i in 1:N){
      # Effect of this subject
      subjEff <- BOLDperc + rnorm(n = 1, mean = 0, sd = bsSD)
      
      # Generate data for each subject in a specific normalization
      normS <- genDat(base = baseJ, BOLDperc = subjEff, 
                      X = X, nscan = nscan, sigmaWN = sigmaWN)
      # Fit glm at subject level
      dataSubj <-
        summary(lm(normS ~ X))$coeff[,'Estimate'] %>% 
        as_data_frame(.) %>%
        mutate(parameter = c('b0', 'b1'),
               subject = i,
               base = baseJ) %>% 
        bind_rows(dataSubj, .)
    }
    # Fit GLM at group level (OLS)
    dataGroup <- 
      dataSubj %>%
      filter(base == baseJ) %>%
      filter(parameter == 'b1') %>%
      lm(value ~ 1, data = .) %>%
      summary(.) %>%
      coef(.) %>%
      # Now select the beta coefficient 
      data.frame(.) %>%
      dplyr::select(Estimate, t.value) %>%
      # Add hedges g
      mutate(hedgeG = NeuRRoStat::hedgeG(t = t.value, N = N, type = 'exact')) %>%
      # Info about the normalization
      mutate(base = baseJ) %>%
      # Simulation
      mutate(sim = k) %>%
      # Add to data frame
      bind_rows(dataGroup, .)
  }
}

# Analysis
dataGroup %>% 
  group_by(base) %>%
  summarise(meanBeta = mean(Estimate),
            meanG = mean(hedgeG))


##
###############
### Generate data: both normalizations
###############
##


# Empty data frame over all simulations
dataGroupB <-
  data.frame() %>% as_tibble()

# For loop over the amount of simulations
for(k in 1:nsim){
  # Empty data frame for subjects
  dataSubj <- data.frame() %>% as_tibble()
  
  # For loop over all subjects: half of them is one normalization
  for(i in 1:N){
    # Switch subjects
    if(i <= N/2) baseJ <- base1
    if(i > N/2) baseJ <- base2
    
    # Effect of this subject
    subjEff <- BOLDperc + rnorm(n = 1, mean = 0, sd = bsSD)
    
    # Generate data for each subject in a specific normalization
    normS <- genDat(base = baseJ, BOLDperc = subjEff, 
                    X = X, nscan = nscan, sigmaWN = sigmaWN)
    # Fit glm at subject level
    dataSubj <-
      summary(lm(normS ~ X))$coeff[,'Estimate'] %>% 
      as_data_frame(.) %>%
      mutate(parameter = c('b0', 'b1'),
             subject = i,
             base = baseJ) %>% 
      bind_rows(dataSubj, .)
  }
  # Loop over both baselines
  for(j in 1:2){
    # Take the corresponding baseline
    if(j == 1) baseJ <- base1
    if(j == 2) baseJ <- base2
    
    # Fit GLM at group level (OLS)
    dataGroupB <- 
    dataSubj %>%
      filter(base == baseJ) %>%
      filter(parameter == 'b1') %>%
      lm(value ~ 1, data = .) %>%
      summary(.) %>%
      coef(.) %>%
      # Now select the beta coefficient 
      data.frame(.) %>%
      dplyr::select(Estimate, Std..Error, t.value) %>%
      # Rename
      rename(Estimate = Estimate, SE = Std..Error, tval = t.value) %>%
      # Add information about the baseline
      mutate(base = baseJ) %>%
      # Scale if baseline == base1
      mutate(Estimate = ifelse(base == base1,
                    Estimate/100,
                    Estimate),
             SE = ifelse(base == base1,
                               SE/100,
                               SE),
             tval = ifelse(base == base1,
                               Estimate/SE,
                               tval)) %>%
      # Add hedges g
      mutate(hedgeG = NeuRRoStat::hedgeG(t = tval, N = N, type = 'exact')) %>%
      # Simulation
      mutate(sim = k) %>%
      # Add to data frame
      bind_rows(dataGroupB, .)
  }
}

# Analysis
dataGroupB %>% 
  group_by(base) %>%
  summarise(meanBeta = mean(Estimate),
            meanG = mean(hedgeG))
































# Quick test
checkSD <- NULL
for(i in 1 :10000){
  test <- X + rnorm(n = nscan, sd = sigmaWN)
  checkSD <- c(checkSD, sd(test))
}
mean(checkSD)
sigmaWN




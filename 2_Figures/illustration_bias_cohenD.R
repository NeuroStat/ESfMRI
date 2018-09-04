####################
#### TITLE:   Illustration of bias in Cohen d
#### Contents:
####
#### Source Files:
#### First Modified: 06/02/2018
#### Notes:
#################

##
###############
### Notes
###############
##

# Small illustration to show the induced bias using Cohen's d and
# both the expressions.


##
###############
### Preparation
###############
##

# Libraries
library(tidyverse)
library(NeuRRoStat)

# Extra functions: first correction factor h
corrH <- function(N){
  num <- gamma((N - 1)/2)
  denom <- (sqrt((N - 1)/2) * gamma((N - 2)/2))
  return(num/denom)
}

##
###############
### Illustration
###############
##

# Create data frame with the number of degrees of freedom and
# the bias using Cohen d (which is just J or h).
CorrFactor <- data.frame('df' = 3:100, 
                         'CorrectionF' = rep(c('J', 'h'), each = length(3:100)),
                         'Value' = c(NeuRRoStat::corrJ(3:100),
                                     corrH(3:100)))
# Using grids
ggplot(CorrFactor, aes(x = df, y = Value)) + 
  geom_line(colour = 'black', size = 1) +
  #geom_line(colour = '#0570b0', size = 1) +
  scale_x_continuous('Degrees of freedom', 
                     minor_breaks = seq(10, 200, by = 5)) +
  scale_y_continuous(expression(delta/E(d))) +
  # ggtitle("Induced bias using Cohen's d",
  #         subtitle = 'True expectation (h) versus approximation (J)') +
  facet_grid(~ CorrectionF) + 
  theme_bw()

# In one plot
ggplot(CorrFactor, aes(x = df, y = Value, group = CorrectionF)) + 
  geom_line(size = 0.5, aes(colour = CorrectionF)) +
  scale_x_continuous('Degrees of freedom', 
                     minor_breaks = seq(10, 200, by = 5)) +
  scale_y_continuous(expression(delta/E(d))) +
  theme_bw()

# Zoom in
CorrFactor %>% 
  filter(df <= 7) %>%
ggplot(., aes(x = df, y = Value, group = CorrectionF)) + 
  geom_line(size = 1, aes(linetype = CorrectionF)) +
  scale_x_continuous('Degrees of freedom', 
                     minor_breaks = seq(10, 200, by = 5)) +
  scale_y_continuous(expression(delta/E(d))) +
  scale_linetype_manual('Correction factor', labels = c('h', 'J'),
                        values = c('solid', 'dashed')) +
  theme_bw() +
  theme(legend.position = 'right')
  
  
# Split window plot
bind_rows(
  data.frame(filter(CorrFactor, CorrectionF == 'h'),
             'window' = 'df = 100', stringsAsFactors = FALSE),
  CorrFactor %>% 
    filter(df <= 7) %>%
    mutate('window' = 'df = 7')
) %>%
  ggplot(., aes(x = df, y = Value, group = CorrectionF)) + 
  geom_line(size = .8, aes(linetype = CorrectionF)) +
  scale_x_continuous('Degrees of freedom', 
                     minor_breaks = seq(10, 200, by = 5)) +
  scale_y_continuous(expression(delta/E(d))) +
  scale_linetype_manual('correction factor', labels = c('h', 'J'),
                        values = c('solid', 'dashed')) +
  facet_wrap(~ window, scales = 'free_x') +
  theme_bw() +
  theme(legend.position = 'bottom')


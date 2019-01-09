####################
#### TITLE:   Plot the different distributions for estimated effects and between-study variability
#### Contents:
####
#### Source Files: //Data/BSVar
#### First Modified: 12/01/2018
#### Notes:
#################


##
###############
### Notes
###############
##

# The RMarkdown file is too heavy to keep on knitting and working on the 
#   different plots.
# Hence I saved the .rds object and will create plots here.
# Then paste them back into the Markdown.



##
###############
### Preparation
###############
##

# Working directory should be source location

# Libraries
library(tidyverse)

# Extra function to create split violin plots
#       source: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <-
  ggproto(
    "GeomSplitViolin",
    GeomViolin,
    draw_group = function(self, data, ..., draw_quantiles = NULL) {
      data <-
        transform(
          data,
          xminv = x - violinwidth * (x - xmin),
          xmaxv = x + violinwidth * (xmax - x)
        )
      grp <- data[1, 'group']
      newdata <-
        plyr::arrange(transform(data, x = if (grp %% 2 == 1)
          xminv
          else
            xmaxv), if (grp %% 2 == 1)
              y
          else-y)
      newdata <-
        rbind(newdata[1,], newdata, newdata[nrow(newdata),], newdata[1,])
      newdata[c(1, nrow(newdata) - 1, nrow(newdata)), 'x'] <-
        round(newdata[1, 'x'])
      if (length(draw_quantiles) > 0 &
          !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                  1))
        quantiles <-
          ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <-
          data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin",
                         grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
      }
      else {
        ggplot2:::ggname("geom_split_violin",
                         GeomPolygon$draw_panel(newdata, ...))
      }
    }
  )



##
###############
### Read in the data
###############
##

estTau <- readRDS('1_Data/2_Processed/estTau.rds')

# Switch to factor
estTau$type <- factor(estTau$type)

##
###############
### Create figures
###############
##


# Begin with the different distributions (weighted average and 
#   between-study heterogeneity) using DL and REML.
paramDL_HE <- estTau %>% 
  filter(type %in% c('DL', 'HE')) %>%
  gather(key = 'Parameter', value = 'Value', -type) %>%
  ggplot(., aes(x = Parameter, y = Value)) +
  geom_violin(aes(fill = type), alpha = 0.5) +
  scale_fill_manual('Between-study heterogeneity estimator: \n', 
                    values = c('#525252', '#d9d9d9')) +
  facet_wrap(type ~ Parameter, scales = 'free', ncol = 2, dir = 'v') +
  coord_flip() +
  scale_x_discrete('Parameter', labels = c('beta' = expression(beta),
                                           'H2'   = expression(H^2),
                                           'I2'   = expression(I^2),
                                           'tau2'   = expression(tau^2))) +
  ggtitle('Observed distributions of various parameters inside ROI') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
paramDL_HE

# Split violin
paramDL_HE_split <- estTau %>% 
  filter(type %in% c('DL', 'HE')) %>%
  gather(key = 'Parameter', value = 'Value', -type) %>%
  ggplot(., aes(x = Parameter, y = Value, fill = type)) +
  geom_split_violin(alpha = 0.4) +
  scale_fill_brewer('', type = 'qual', palette = 2) +
  facet_wrap( ~ Parameter, scales = 'free', ncol = 2, dir = 'v') +
  coord_flip() +
  scale_x_discrete('Parameter', labels = c('beta' = expression(hat(beta)),
                                           'H2'   = expression(H^2),
                                           'I2'   = expression(I^2),
                                           'tau2'   = expression(hat(tau)^2))) +
  ggtitle('Empirical distributions of various parameters within ROI') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
paramDL_HE_split


# Now plot between-study heterogeneity using all estimators
allEst <- estTau %>% 
  dplyr::select(tau2, type) %>%
  ggplot(., aes(x = type, y = tau2)) +
  geom_violin(aes(fill = type), alpha = 0.5, trim = TRUE) +
  scale_fill_brewer('', type = 'qual', palette = 2) +
  facet_wrap(~ type, scales = 'free_y', ncol = 1) +
  coord_flip() +
  scale_x_discrete('Estimator') +
  scale_y_continuous(expression(hat(tau)^2)) +
  guides(fill = FALSE)  + 
  ggtitle('Empirical distributions of between-study heterogeneity') +
  labs(caption = 'Over voxels within ROI') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
allEst


# Plot estimated effect sizes using all estimators
allBeta_est <- estTau %>% 
  dplyr::select(beta, type) %>%
  ggplot(., aes(x = type, y = beta)) +
  geom_violin(aes(fill = type), alpha = 0.5, trim = TRUE,
              draw_quantiles = c(.95),
              linetype = 'solid') +
  scale_fill_brewer('', type = 'qual', palette = 2) +
  facet_wrap(~ type, scales = 'free_y', ncol = 1) +
  coord_flip() +
  scale_x_discrete('Estimator') +
  scale_y_continuous(expression(hat(beta)^2)) +
  guides(fill = FALSE)  + 
  ggtitle('Empirical distributions of estimated effect size') +
  labs(caption = 'Over voxels within ROI') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
allBeta_est


##
###############
### Save figures
###############
##


paramDL_HE
paramDL_HE_split
allEst
allBeta_est

# Save various parameters using either DL or HE
ggsave(filename = paste0(getwd(), '/2_Figures/paramDL_HE.png'),
       plot = paramDL_HE,
       width = 20, height = 24, units = 'cm', scale = 1)

# Save various parameters using either DL or HE: split violin
ggsave(filename = paste0(getwd(), '/2_Figures/paramDL_HE_split.png'),
       plot = paramDL_HE_split,
       width = 22, height = 22, units = 'cm', scale = 1)

# Save all beta distributions
ggsave(filename = paste0(getwd(), '/2_Figures/allBeta_est.png'),
       plot = allBeta_est,
       width = 22, height = 26, units = 'cm', scale = 1)

# Save all between-study heterogeneity estimators
ggsave(filename = paste0(getwd(), '/2_Figures/allEst.png'),
       plot = allEst,
       width = 22, height = 26, units = 'cm', scale = 1)










###########################################################################
# WORKSPACE PREP
###########################################################################

rm(list=ls())

# Set trees/burns/draws for tsbart and bart.
ntree=200; nburn=100; nsim=250

#-------------------------------------------------------------------------
# Install tsbart and fastbart packages.
#-------------------------------------------------------------------------

library(tsbart)
library(dbarts)
library(splines)
library(caret)       # For test/train split
library(data.table)  # For fast csv read/write.
library(gridExtra)   # For plotting multi-pane ggplots.
library(tidyverse)   # For tidy data/plotting.
library(mgcv)        # For penalized splines
library(magrittr)
library(cowplot)
library(ggthemes)
library(viridis)

# Set output directories.
out.fig = './output-figures/'
out.csv = './output-files/'

# Load custom functions.
source('./code/helper-functions/LBW-functions-modelfitutils.R')
source('./code/helper-functions/LBW-functions-cvutils.R')
source('./code/helper-functions/LBW-functions-testtrain-split.R')
source('./code/helper-functions/LBW-functions-binll.R')
source('./code/helper-functions/ggtheme-publication.R')

# ggplot settings
theme_set(theme_bw(base_size=16, base_family='Helvetica'))

###########################################################################
# DATA LOAD AND PREP
###########################################################################

#-------------------------------------------------------------------
# 1. Test/train data from case-control sampled obstetrics dataset.

# Read rds file and expand data frame for correct conditional probabilities.
data = as.data.frame(readRDS(paste0(getwd(),'/data/mainland_sample.rds'))) %>%
  mutate(LBW = as.numeric(LBW) -1)
data = upsample(data, n=2000, sb_pct = .5)
data = tsbart::survPrep(data, 'combgest', 'LBW') 

#-------------------------------------------------------------------
# 2. Set the column names of covariates to include in the models.
xcols = c("MRACE6", "mager", "RF_GDIAB", "RF_GHYPE", "preterm",
          "previs", "wtgain", "CIG_REC", "RF_PPTERM", "SEX")

#-------------------------------------------------------------------
# 3. Test-train split.
set.seed(64)
train_idx <- createDataPartition(data$LBW, p = 0.8, list = FALSE)
train <- data[train_idx, ]
test <- data[-train_idx, ]
test$testtrain = 'test'
train$testtrain = 'train'

#-------------------------------------------------------------------
# 3. Read in new patient panel.
newpts = readRDS(paste0(getwd(), '/data/territory_sample.rds')) %>%
  mutate(LBW = as.numeric(as.character(LBW))) 
newpts$id = newpts$id + 5000000  #Panel IDs start at 5000001.
ptpanel = tsbart::makePredGrid(newpts, 'combgest', sort(unique(data$combgest)))  

ptpanel$testtrain = 'panel'

# Add patient panel to end of training data.
#    (Want models to predict on both the test data and the new patient panel,
#    for plotting and OOS log-loss.)

test = dplyr::bind_rows(ptpanel[c(xcols, 'id', 'combgest','LBW','testtrain')],
                        test[c(xcols, 'id', 'combgest','LBW','testtrain')]) %>% 
  filter(RF_GDIAB != "Unknown")

#test$testtrain = c(rep('panel', nrow(ptpanel)), rep('test',nrow(test) - nrow(ptpanel)))


###########################################################################
# Fit models.
###########################################################################

# tsbart with previously selected tuned hyperparameter of 0.6
cv_tsb_panel = tsbart_fit_util(train, test, xcols, ec=0.6, ntree=ntree, nburn=nburn, nsim=nsim)

# tsbart with default hyperparameter 1
#cv_tsb_default_panel = tsbart_fit_util(train, test, xcols, ec=1, ntree=ntree, nburn=nburn, nsim=nsim)

# vanilla bart
cv_vb_panel = bart_fit_util(train, test, xcols=c(xcols,'combgest'), ntree=ntree, nburn=nburn, nsim=nsim)

# Fit splines model with linear interaction.
cv_sp1_panel = spline_fit_util(train, test)

# Fit splines model with basis interaction.
cv_sp2_panel = spline_fit_util(train, test, spline_interactions=TRUE)

# Fit p-splines model.
cv_pensp_panel = penspline_fit_util(train, test)

# Fit random forest model.
cv_rf_panel = rf_fit_util(train, test)

# Extract fits.
test_panel = extractFits(test, 
                   cv_tsb = cv_tsb_panel, 
                   cv_vb = cv_vb_panel, 
                   cv_sp1 = cv_sp1_panel, 
                   cv_sp2 = cv_sp2_panel, 
                   cv_pensp = cv_pensp_panel, 
                   rf = cv_rf_panel, 
                   cv_tsb_default=NULL, 
                   adjust=T)
ptpanel = test_panel %>% filter(testtrain=='panel')


###########################################################################
# Plotting
###########################################################################

y_limits = c(0,10)
x_limits = c(25, 42)
ids_newplt = c(5000013, 5000119, 5000173, 5000177, 5000182)

bkg_linetype = 2

row1 = panplt_tsb_combos = ggplot((ptpanel %>% filter(id %in% ids_newplt) %>% arrange(id)), 
                                  aes(x=combgest, y=phat_oos_tsb), colour='black') +
   geom_line(aes(y=phat_oos_vb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp1), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp2), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_pensp), colour='grey60', linetype=bkg_linetype ) +
   geom_ribbon(aes(x=combgest, ymin=phat_oos_tsb_lb, ymax=phat_oos_tsb_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line(size=.75,colour='blue') +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   coord_cartesian(ylim=y_limits,  xlim = x_limits) +
   labs(x='',y='tsBART') +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines")) +
      cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85') +
   cowplot::panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE) +
   theme(plot.margin = unit(c(0, .2, 0, .05), "cm"))

row2 = ggplot((ptpanel %>% filter(id %in% ids_newplt) %>% arrange(id)),
              aes(x=combgest, y=phat_oos_vb), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp1), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp2), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_pensp), colour='grey60', linetype=bkg_linetype ) +
   geom_ribbon(aes(x=combgest, ymin=phat_oos_vb_lb, ymax=phat_oos_vb_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line(size=.75,colour='blue') +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='BART') +
   coord_cartesian(ylim=y_limits, xlim = x_limits) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines")) +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85') +
   cowplot::panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE)+
   theme(plot.margin = unit(c(0, .2, 0, .05), "cm"))

row3 = ggplot((ptpanel %>% filter(id %in% ids_newplt) %>% arrange(id)),
              aes(x=combgest, y=phat_oos_sp1), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_vb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp2), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_pensp), colour='grey60', linetype=bkg_linetype ) +
   geom_ribbon(aes(x=combgest, ymin=phat_oos_sp1_lb, ymax=phat_oos_sp1_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line(size=.75,colour='blue') +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='Splines') +
   coord_cartesian(ylim=y_limits, xlim = x_limits) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines"))  +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85') +
   cowplot::panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE)+
   theme(plot.margin = unit(c(0, .2, 0, .05), "cm"))

# row4 = ggplot((ptpanel %>% filter(id %in% ids_newplt)), aes(x=combgest, y=phat_oos_sp2), colour='black') +
#    geom_line(aes(y=phat_oos_tsb), colour='grey60', linetype=bkg_linetype ) +
#    geom_line(aes(y=phat_oos_vb), colour='grey60', linetype=bkg_linetype ) +
#    geom_line(aes(y=phat_oos_sp1), colour='grey60', linetype=bkg_linetype ) +
#    geom_line(aes(y=phat_oos_pensp), colour='grey60', linetype=bkg_linetype ) +
#    geom_ribbon(aes(x=combgest, ymin=phat_oos_sp2_lb, ymax=phat_oos_sp2_ub), alpha=0.35, fill='dodgerblue4') +
#    geom_line(size=.75,colour='blue') +
#    facet_wrap(~label, ncol=length(ids_newplt)) +
#    labs(x = '',y='Splines 2') +
#    coord_cartesian(ylim=y_limits) +
#    scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines"))  +
#    theme(
#       strip.background = element_blank(),
#       strip.text.x = element_blank(),
#       axis.text.x=element_blank()
#    ) +
#    cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85') +
#    cowplot::panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE)+
#    theme(plot.margin = unit(c(0, .2, 0, .05), "cm"))

row5 = ggplot((ptpanel %>% filter(id %in% ids_newplt) %>% arrange(id)),
              aes(x=combgest, y=phat_oos_pensp), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_vb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp1), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp2), colour='grey60', linetype=bkg_linetype ) +
   geom_ribbon(aes(x=combgest, ymin=phat_oos_pensp_lb, ymax=phat_oos_pensp_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line(size=.75,colour='blue') +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='P-splines') +
   coord_cartesian(ylim=y_limits, xlim = x_limits) +
   scale_x_continuous(expand=c(0,0)) + theme(panel.spacing.x=unit(1.25,"lines"))  +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x=element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85') +
   cowplot::panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE)+
   theme(plot.margin = unit(c(0, .2, 0, .05), "cm"))

row6 = ggplot((ptpanel %>% filter(id %in% ids_newplt) %>% arrange(id)),
              aes(x=combgest, y=phat_oos_rf), colour='black') +
   geom_line(aes(y=phat_oos_tsb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_vb), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp1), colour='grey60', linetype=bkg_linetype ) +
   geom_line(aes(y=phat_oos_sp2), colour='grey60', linetype=bkg_linetype ) +
   geom_ribbon(aes(x=combgest, ymin=phat_oos_rf_lb, ymax=phat_oos_rf_ub), alpha=0.35, fill='dodgerblue4') +
   geom_line(size=.75,colour='blue') +
   facet_wrap(~label, ncol=length(ids_newplt)) +
   labs(x = '',y='Random Forest') +
   coord_cartesian(ylim=y_limits, xlim = x_limits) +
  scale_x_continuous(
    breaks = c(25, 30, 35, 40),  # Explicitly set breaks at 25, 30, 35
    labels = c("25", "30", "35", "40"),  # Optional: Ensures clean labels
    expand = c(0, 0)) + 
      theme(panel.spacing.x=unit(1.25,"lines"))  +
   theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
   ) +
   cowplot::background_grid(major = "xy", minor = "none", colour.major = 'grey85') +
   cowplot::panel_border(colour = "black", size = 0.5, linetype = 1, remove = FALSE) +
   theme(plot.margin = unit(c(0, .2, 0, .05), "cm"))

library(grid)
grid.newpage()
all = grid.arrange(row1, row2, row3, row5, row6,
             ncol=1,
             heights = c(1, 1, 1, 1, 1.25),
             bottom='Gestational age (Wks)',
               left='Risk of LBW per 100 remaining pregnancies',
             top='            Patient 1                 Patient 2                 Patient 3                Patient 4                Patient 5')

ggsave(paste0(out.fig,'figure-panel-prediction.png'), all,
       width=8, height=8, units='in', dpi=300, limitsize=TRUE)


save.image("./`Final Paper`.RData")

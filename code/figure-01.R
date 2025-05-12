# For creating Figure 1.
# Note: requires tsbart package to reproduce regular BART code.
# tsbart package can be obtained by request from Jared S. Murray.

rm(list=ls())

library(Rcpp)
library(RcppArmadillo)

#===================================================================
# Workspace prep
#===================================================================

# Load other libraries.
library(tsbart)
library(dbarts)
library(mosaic)
library(tidyverse)
library(gridExtra)


# ggplot settings
source('./code/helper-functions/ggtheme-publication.R')
theme_set(theme_bw(base_size=16, base_family='Helvetica'))

# Output directories.
out.fig = paste0(getwd(),'/output-figures/')


##################################################################################
### Loop through each simulation file, if an individual name is not specified.
##################################################################################

# Optimize expected number of crossings over user-defined grid.
exp_cross = NULL

#===================================================================
# Read data.
#===================================================================

#-------------------------------------------------------------------
# 1. Test/train data.

# Read csv file and expand data frame for correct conditional probabilities.
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

# Extract values from data frame.
yob = train$LBW;       y_pred = test$LBW;      y= ifelse(yob==0, -1.96, 1.96)        # Vector of responses.
ti = train$combgest;     ti_pred = test$combgest      # Time points for each obs.
#fx = train$fx;     fx_pred = test$fx      # True underlying function, without epsilon noise.

xx = train[xcols]; x_pred = test[xcols]     # Matrix of predicted covariates.


#############################################################################################
###   1. tsBART.
#############################################################################################

#=====================================================================
#=== Evaluate optimal expected number of crossings.
#=====================================================================

# Calibrate BART's error variance a la CGM 2010 (method 2)
df = data.frame(xx, 't'=train$combgest, 'y'=train$LBW)
lmf = lm(y~., df)
sighat = sigma(lmf) 

# Hyperparameters
nu = 3
sigq = .9
qchi = qchisq(1.0-sigq,nu)
lambda = (sighat*sighat*qchi)/nu

# Evaluate optimal number of crossings.
ecross_candidates = seq(0.1, 1,by=.1)
ecrossTune = tuneEcross(ecross_candidates,
                        y=y, tgt=ti, x=xx, 
                        nburn=200, nsim=1000, ntree=200,
                        sigq=sigq, nu=nu,
                        base_tree=.95, power_tree=2,
                        probit=TRUE, yobs=yob)


# Set expected number of crossings.
exp_cross = ecrossTune$ecross_opt
waic_plt = ecrossTune$waic_plot + theme_Publication()


# Output figure for appendix.
ggsave(paste0(out.fig,"figure-tuned-crossing.png"),
       waic_plt, height=5, width=7, dpi=300)


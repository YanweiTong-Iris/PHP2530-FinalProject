###########################################################################
# WORKSPACE PREP
###########################################################################

rm(list=ls())

# Set trees/burns/draws for tsbart and bart.
ntree=200; nburn=1000; nsim=5000

#-------------------------------------------------------------------------
# Workspace setup.
#-------------------------------------------------------------------------

library(tsbart)
library(dbarts)
library(splines)

# Load other libraries.
library(data.table)  # For fast csv read/write.
library(gridExtra)   # For plotting multi-pane ggplots.
library(tidyverse)   # For tidy data/plotting.
library(mgcv)        # For penalized splines
library(viridis)
library(ggthemes)
library(caret)

# ggplot settings
theme_set(theme_bw(base_size=16, base_family='Helvetica'))

# Load helper functions.
source('./code/helper-functions/LBW-functions-modelfitutils.R')
source('./code/helper-functions/LBW-functions-cvutils.R')
source('./code/helper-functions/LBW-functions-testtrain-split.R')
source('./code/helper-functions/LBW-functions-binll.R')
source('./code/helper-functions/ggtheme-publication.R')

###########################################################################
# DATA LOAD AND PREP
###########################################################################

#-------------------------------------------------------------------
# 1. Test/train data from case-control sampled obstetrics dataset.

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
train_idx <- createDataPartition(data$preterm, p = 0.8, list = FALSE)
train <- data[train_idx, ]
test <- data[-train_idx, ]
test$testtrain = 'test'
train$testtrain = 'train'

###########################################################################
# Fit models.
###########################################################################

# tsbart with previously tuned hyperparameter of .1.
cv_tsb = tsbart_fit_util(train, test, xcols, ec=5, ntree=ntree, nburn=nburn, nsim=nsim, init_upper=.01, init_lower=qnorm(.05))

# tsbart with default hyperparameter 1
#cv_tsb_default = tsbart_fit_util(train, test, xcols, ec=1, ntree=ntree, nburn=nburn, nsim=nsim, init_upper=.01, init_lower=qnorm(.05))

# vanilla bart
cv_vb = bart_fit_util(train, test, xcols=c(xcols,'combgest'), ntree=ntree, nburn=nburn, nsim=nsim, init_upper=.01, init_lower=qnorm(.05))

# Fit splines model with linear interaction.
cv_sp1 = spline_fit_util(train, test)

# Fit splines model with basis interaction.
cv_sp2 = spline_fit_util(train, test, spline_interactions=FALSE)

# Fit p-splines model.
cv_pensp = penspline_fit_util(train, test)

# Fit random forest model.
cv_rf = rf_fit_util(train, test)


# Add out of sample fit info to dataframe.
test_oos = extractFits(test, 
                   cv_tsb = cv_tsb, 
                   cv_vb = cv_vb, 
                   cv_sp1 = cv_sp1, 
                   cv_sp2 = cv_sp2, 
                   cv_pensp = cv_pensp, 
                   rf = cv_rf, 
                   cv_tsb_default=NULL, 
                   adjust=T)

###########################################################################
# Binomial log-likelihoods (per person) & Weekly misclassification error.
###########################################################################
library(pROC)
roc(response = test_oos$LBW, predictor = test_oos$phat_oos_tsb)
roc(response = test_oos$LBW, predictor = test_oos$phat_oos_vb)
roc(response = test_oos$LBW, predictor = test_oos$phat_oos_sp1)
roc(response = test_oos$LBW, predictor = test_oos$phat_oos_pensp)

roc_obj <- roc(response = test_oos$LBW, predictor = test_oos$phat_oos_tsb)
thrsh = coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
thrsh = thrsh[1]

binll_oos = binll_per_person(test_oos, method='tsb', thresh=thrsh)

for(meth in c('vb','sp1','pensp')){
   binll_oos = rbind.data.frame(binll_oos, binll_per_person(test_oos, meth, thresh=thrsh))
}

# Display results.
binll_oos %>% filter(combgest =='all') %>% dplyr::select(method,logl)
binll_oos %>% filter(combgest !='all')



write.csv(binll_oos %>% 
            filter(combgest =='all') %>% 
            dplyr::select(method,logl), './output-files/table-oos-logloss.csv', row.names=F)


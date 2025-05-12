# Helper functions for fitting each model, for binary response.

########################################################################
# tsbart fit utility.
########################################################################

tsbart_fit_util = function(train, test, xcols, yvar='LBW', tvar='combgest', 
                           ec=1, ntree=200, nburn=200, nsim=2500, padjust=50, 
                           init_upper=1.96, init_lower=-1.96){   
   
   require(tsbart)
   require(msm)

   #===================================================================
   # Extract data for model inputs.
   #===================================================================
   
   phat = mean(train[,paste(yvar)])
   offset = qnorm(phat) 
   
   #-------------------------------------------------------------------
   # Model matrices and cutpoints for tsbart. (Does not include time.)
   #-------------------------------------------------------------------
   
   # With time included.
   xx = tsbart::makeModelMatrix(train[xcols])
   xpred = tsbart::makeModelMatrix(test[xcols])
   cutpoints = tsbart::makeCutpoints(xx)
   
   #-------------------------------------------------------------------
   # Response, latent response, and time points.
   #-------------------------------------------------------------------

   # Training set.
   yobs = unlist(train[paste(yvar)])
   y = ifelse(yobs==1,init_upper, init_lower) 
   
   # Test set.
   yobs_pred = unlist(test[paste(yvar)])
   ypred =  ifelse(yobs==1,init_upper, init_lower) 
   
   # Time points.
   t = unlist(train[paste(tvar)])
   tpred = unlist(test[paste(tvar)])
   
   #===================================================================
   # tsbart setup for model fitting.  (Hyperparameters, etc.)
   #===================================================================
   
   # Calibrate BART's error variance a la CGM 2010 (method 2)
   df = data.frame(xx, t, y)
   lmf = lm(y~.,df)
   sighat = sigma(lmf) 
   
   # Hyperparameters
   nu = 3
   sigq = .9
   qchi = qchisq(1.0-sigq,nu)
   lambda = (sighat*sighat*qchi)/nu
   
   phat = mean(unlist(yobs))
   offset = qnorm(phat)
   
   #===================================================================
   # Tune expected number of crossings.
   #===================================================================
   exp_cross = ec

   if(is.null(exp_cross)){
      ecross_candidates = c(.1,.5,1,2.5,5)

      # Evaluate optimal number of crossings.
      ecrossTune = tuneEcross(ecross_candidates, 
                              y=y, tgt=tgt, x=xx, 
                              nburn=250, nsim=200, ntree=200,
                              yobs=yobs, probit=T)

      # Set expected number of crossings.
      exp_cross = ecrossTune$ecross_opt
   }

   print(exp_cross)
   
   #===================================================================
   # Fit tsbart model.
   #===================================================================

   fit = tsbart(
     y = y,
     yobs = yobs,
     tgt =  t,
     tpred =  tpred,
     x = xx,
     xpred = xpred,
     nburn = nburn,
     nsim = nsim,
     ntree = ntree,
     ecross = exp_cross,
     probit = TRUE
   )
   
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   phat_mcmc_oos = pnorm(fit$mcmcdraws_oos)

   phat_oos = apply(phat_mcmc_oos, 2, mean)
   phat_oos_lb =  apply(phat_mcmc_oos, 2, function(x) quantile(x, 0.025))
   phat_oos_ub =  apply(phat_mcmc_oos, 2, function(x) quantile(x, 0.975))
   
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   #===================================================================
   # Results
   #===================================================================
   
   #-------------------------------------------------------------------
   # Return output list.
   #-------------------------------------------------------------------
   return(list('fit' = fit,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos_lb,
               'phat_oos_ub' = phat_oos_ub,
               'binll_oos' = binll_oos,
               't'=t,
               'tpred'=tpred,
               'exp_cross' = exp_cross))
}


########################################################################
# Vanilla BART fit utility.
########################################################################

bart_fit_util = function(train, test, xcols, yvar='LBW', tvar='combgest', ntree=200, nburn=200, nsim=2500,
                         init_upper=1.96, init_lower=-1.96){   
   
   # Note: xcols should include the time variable here.
   require(dbarts)
   
   #===================================================================
   # Extract data for model inputs.
   #===================================================================
   
   phat = mean(train[,paste(yvar)])
   offset = qnorm(phat)
   
   #-------------------------------------------------------------------
   # Model matrices and cutpoints for tsbart.  With time included.  This is for vanilla BART.
   #-------------------------------------------------------------------
   xx = tsbart::makeModelMatrix(train[xcols])
   xpred = tsbart::makeModelMatrix(test[xcols])
   cutpoints = tsbart::makeCutpoints(xx)
   
   #-------------------------------------------------------------------
   # Response, latent response, and time points.
   #-------------------------------------------------------------------
   
   # Training set.
   yobs = unlist(train[paste(yvar)])
   y = ifelse(yobs==1,init_upper, init_lower) 
   
   # Test set.
   yobs_pred = unlist(test[paste(yvar)])
   ypred =  ifelse(yobs==1,init_upper, init_lower) 
   
   # Time points.
   t = unlist(train[paste(tvar)])
   tpred = unlist(test[paste(tvar)])
   
   #-------------------------------------------------------------------
   # BART setup for model fitting.  (Hyperparameters, etc.)
   #-------------------------------------------------------------------
   
   # Calibrate BART's error variance a la CGM 2010 (method 2)
   df = data.frame(xx, t, y)
   lmf = lm(y~.,df)
   sighat = sigma(lmf)
   
   # Hyperparameters
   nu = 3
   sigq = .9
   qchi = qchisq(1.0-sigq,nu)
   lambda = (sighat*sighat*qchi)/nu
   
   phat = mean(unlist(yobs))
   offset = qnorm(phat)
   
   #-------------------------------------------------------------------
   # Fit BART model.
   #-------------------------------------------------------------------
   # fit = fastbart::bartRcppClean(y_ = y, x_ = t(xx), # obsvervations must be in *columns*
   #                               xpred_ = t(xpred), 
   #                               xinfo_list = cutpoints,
   #                               nburn, nsim, ntree,
   #                               lambda, nu, kfac=2,
   #                               paste0(getwd(),"/vanillabarttrees.txt"),
   #                               RJ=FALSE)
   
   fit = bart(x.train=train[xcols], 
              y.train=train$LBW, 
              x.test=test[xcols],
              ntree=ntree,
              ndpost=nsim,
              nskip=nburn, 
              binaryOffset=offset)
   
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   phat_mcmc_oos = pnorm(fit$yhat.test)
   
   phat_oos = apply(phat_mcmc_oos, 2, mean)
   phat_oos_lb =  apply(phat_mcmc_oos, 2, function(x) quantile(x, 0.025))
   phat_oos_ub =  apply(phat_mcmc_oos, 2, function(x) quantile(x, 0.975))
   
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   #===================================================================
   # Results
   #===================================================================
   
   return(list('fit' = fit,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos_lb,
               'phat_oos_ub' = phat_oos_ub,
               'binll_oos' = binll_oos,
               't'=t,
               'tpred'=tpred))
}


########################################################################
# Spline fit utility (fits preliminary spline model).
########################################################################
spline_fit_util = function(train, test, spline_interactions=F){
   
   #-------------------------------------------------------------------
   # Response, latent response, and time points.
   #-------------------------------------------------------------------
   yobs = train$LBW
   yobs_pred = test$LBW
   
   #-----------------------------------------------------------------------
   # First, collapse data for optimal computing speed.
   #-----------------------------------------------------------------------
   
   # Columns to collapse on.
   coll = c("combgest","MRACE6", "mager", "RF_GDIAB", "RF_GHYPE", 
            "previs", "wtgain", "CIG_REC", "RF_PPTERM", "SEX")
   
   # Collapse train set
   trainc = train[c('LBW',coll)] %>%
     group_by(combgest, MRACE6, mager, RF_GDIAB, RF_GHYPE, 
              previs, wtgain, CIG_REC, RF_PPTERM, SEX) %>%
     summarize(LBW_yes = sum(LBW==1), LBW_total = length(LBW), n = n())
   
   trainc$hr = ifelse(trainc$RF_GDIAB==0 & trainc$RF_GHYPE==0 &  trainc$RF_PPTERM==0, 0, 1)
   
   #-----------------------------------------------------------------------
   # Fit spline model.
   #-----------------------------------------------------------------------
   
   require(splines)
   
   # Set degrees of freedom and degree for all spline models.
   degr = 3
   degf = 7 
   
   fitmod=NULL
   
   # Fit training set.
   if(spline_interactions==FALSE){
      fitmod = glm(cbind(LBW_yes, LBW_total) ~ bs(combgest, degree=degr, df=degf) +
                     MRACE6 + mager + RF_GDIAB + RF_GHYPE + previs+ wtgain + CIG_REC + RF_PPTERM + SEX,
                   data = trainc, family='binomial')
      
   } 
   if(spline_interactions==TRUE){
      fitmod = glm(cbind(LBW_yes, LBW_total) ~ 
                      bs(combgest, degree=degr, df=degf):MRACE6 +
                      bs(combgest, degree=degr, df=degf):mager +
                      bs(combgest, degree=degr, df=degf):RF_GDIAB +
                      bs(combgest, degree=degr, df=degf):RF_GHYPE +
                      bs(combgest, degree=degr, df=degf):RF_PPTERM +
                      bs(combgest, degree=degr, df=degf):SEX +
                     bs(combgest, degree=degr, df=degf):CIG_REC +
                      bs(combgest, degree=degr, df=degf):previs + 
                      bs(combgest, degree=degr, df=degf):wtgain,
                   data = trainc, family='binomial')
   }
   
   # Predict on training set.
   predict_is = predict(fitmod, newdata=train, type='response', se=T)
   phat_is = as.vector(as.numeric(predict_is$fit)) 
   se_is = as.vector(as.numeric(predict_is$se.fit)) 
   
   # Predict on test set.
   predict_oos = predict(fitmod, newdata=test, type='response', se=T)
   phat_oos = as.vector(as.numeric(predict_oos$fit))
   se_oos = as.vector(as.numeric(predict_oos$se.fit))  
   
   fit = list('phat_is' = phat_is, 'se_is' = se_is, 'phat_oos' = phat_oos, 'se_oos' = se_oos)
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   return(list('fit' = fit,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos - 1.96 * se_oos,
               'phat_oos_ub' = phat_oos + 1.96 * se_oos,
               'binll_oos' = binll_oos))
}

########################################################################
# Spline fit utility (fits preliminary spline model).
########################################################################
spline_fit_util_sim = function(train, test, spline_interactions=F){
   
   #-------------------------------------------------------------------
   # Response, latent response, and time points.
   #-------------------------------------------------------------------
   yobs = train$LBW
   yobs_pred = test$LBW
   
   #-----------------------------------------------------------------------
   # Fit spline model.
   #-----------------------------------------------------------------------
   
   require(splines)
   
   # Set degrees of freedom and degree for all spline models.
   degr = 3
   degf = 7 
   
   fitmod=NULL
   
   # Fit training set.
   if(spline_interactions==FALSE){
      fitmod = glm(y ~ bs(combgest, degree=degr, df=degf) +
                     MRACE6 + mager + RF_GDIAB + RF_GHYPE + previs+ wtgain + CIG_REC + RF_PPTERM + SEX,
                   data = train, family='binomial')
      
   } 
   if(spline_interactions==TRUE){
     fitmod = glm(cbind(LBW_yes, LBW_total) ~ 
                    bs(combgest, degree=degr, df=degf):MRACE6 +
                    bs(combgest, degree=degr, df=degf):mager +
                    bs(combgest, degree=degr, df=degf):RF_GDIAB +
                    bs(combgest, degree=degr, df=degf):RF_GHYPE +
                    bs(combgest, degree=degr, df=degf):RF_PPTERM +
                    bs(combgest, degree=degr, df=degf):SEX +
                    bs(combgest, degree=degr, df=degf):CIG_REC +
                    bs(combgest, degree=degr, df=degf):previs + 
                    bs(combgest, degree=degr, df=degf):wtgain,
                  data = train, family='binomial')
   }
   
   # Predict on training set.
   predict_is = predict(fitmod, newdata=train, type='response', se=T)
   phat_is = as.vector(as.numeric(predict_is$fit)) 
   se_is = as.vector(as.numeric(predict_is$se.fit)) 
   
   # Predict on test set.
   predict_oos = predict(fitmod, newdata=test, type='response', se=T)
   phat_oos = as.vector(as.numeric(predict_oos$fit))
   se_oos = as.vector(as.numeric(predict_oos$se.fit))  
   
   fit = list('phat_is' = phat_is, 'se_is' = se_is, 'phat_oos' = phat_oos, 'se_oos' = se_oos)
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   return(list('fit' = fit,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos - 1.96 * se_oos,
               'phat_oos_ub' = phat_oos + 1.96 * se_oos,
               'binll_oos' = binll_oos))
}

########################################################################
# Penalized spline fit utility.
########################################################################
penspline_fit_util = function(train, test){
   
   #-------------------------------------------------------------------
   # Response, latent response, and time points.
   #-------------------------------------------------------------------
   yobs = train$LBW
   yobs_pred = test$LBW
   
   #-----------------------------------------------------------------------
   # First, collapse data for optimal computing speed.
   #-----------------------------------------------------------------------
   
   # Columns to collapse on.
   coll = c("combgest","MRACE6", "mager", "RF_GDIAB", "RF_GHYPE", 
            "previs", "wtgain", "CIG_REC", "RF_PPTERM", "SEX")

   
   # Collapse train set
   trainc = train[c('LBW',coll)] %>%
     group_by(combgest, MRACE6, mager, RF_GDIAB, RF_GHYPE, 
              previs, wtgain, CIG_REC, RF_PPTERM, SEX) %>%
      summarize(LBW_yes = sum(LBW==1), LBW_total = length(LBW), n = n())
   
   trainc$hr = ifelse(trainc$RF_GDIAB==0 & trainc$RF_GHYPE==0 &  trainc$RF_PPTERM==0, 0, 1)
   
   #-----------------------------------------------------------------------
   # Fit spline model.
   #-----------------------------------------------------------------------
   
   require(mgcv)
   
   # These are cubic p-splines, with 9 basis elements, with a second-order penalty.
   # See: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.terms.html
   fitmod = gam(cbind(LBW_yes, LBW_total) ~ s(combgest, bs='ps', m=c(2,2)) +
                  MRACE6 + mager + RF_GDIAB + RF_GHYPE + previs+ wtgain + CIG_REC + RF_PPTERM + SEX,
                data = trainc, family='binomial')
   
   # Predict on training set.
   predict_is = predict(fitmod, newdata=train, type='response', se=T)
   phat_is = as.vector(as.numeric(predict_is$fit)) 
   se_is = as.vector(as.numeric(predict_is$se.fit))
   
   # Predict on test set.
   predict_oos = predict(fitmod, newdata=test, type='response', se=T)
   phat_oos = as.vector(as.numeric(predict_oos$fit))
   se_oos = as.vector(as.numeric(predict_oos$se.fit))  
   
   fit = list('phat_is' = phat_is, 'se_is' = se_is, 'phat_oos' = phat_oos, 'se_oos' = se_oos)
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   return(list('fit' = fit,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos - 1.96 * se_oos,
               'phat_oos_ub' = phat_oos + 1.96 * se_oos,
               'binll_oos' = binll_oos))
}

penspline_fit_util_sim = function(train, test){
   
   #-------------------------------------------------------------------
   # Response, latent response, and time points.
   #-------------------------------------------------------------------
   yobs = train$LBW
   yobs_pred = test$LBW
   
   #-----------------------------------------------------------------------
   # Fit spline model.
   #-----------------------------------------------------------------------
   
   require(mgcv)
   
   # These are cubic p-splines, with 9 basis elements, with a second-order penalty.
   # See: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.terms.html
   fitmod = gam(y ~ s(combgest, bs='ps', m=c(2,2)) +
                  MRACE6 + mager + RF_GDIAB + RF_GHYPE + previs+ wtgain + CIG_REC + RF_PPTERM + SEX,
                data = train, family='binomial')
   
   # Predict on training set.
   predict_is = predict(fitmod, newdata=train, type='response', se=T)
   phat_is = as.vector(as.numeric(predict_is$fit)) 
   se_is = as.vector(as.numeric(predict_is$se.fit))
   
   # Predict on test set.
   predict_oos = predict(fitmod, newdata=test, type='response', se=T)
   phat_oos = as.vector(as.numeric(predict_oos$fit))
   se_oos = as.vector(as.numeric(predict_oos$se.fit))  
   
   fit = list('phat_is' = phat_is, 'se_is' = se_is, 'phat_oos' = phat_oos, 'se_oos' = se_oos)
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   return(list('fit' = fit,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos - 1.96 * se_oos,
               'phat_oos_ub' = phat_oos + 1.96 * se_oos,
               'binll_oos' = binll_oos))
}



########################################################################
# Penalized spline fit utility.
########################################################################
penspline_fit_util = function(train, test){
   
   #-------------------------------------------------------------------
   # Response, latent response, and time points.
   #-------------------------------------------------------------------
   yobs = train$LBW
   yobs_pred = test$LBW
   
   #-----------------------------------------------------------------------
   # First, collapse data for optimal computing speed.
   #-----------------------------------------------------------------------
   
   # Columns to collapse on.
   coll = c("combgest","MRACE6", "mager", "RF_GDIAB", "RF_GHYPE", 
            "previs", "wtgain", "CIG_REC", "RF_PPTERM", "SEX")
   
   # Collapse train set
   trainc = train[c('LBW',coll)] %>%
     group_by(combgest, MRACE6, mager, RF_GDIAB, RF_GHYPE, 
              previs, wtgain, CIG_REC, RF_PPTERM, SEX) %>%
     summarize(LBW_yes = sum(LBW==1), LBW_total = length(LBW), n = n())
   
   trainc$hr = ifelse(trainc$RF_GDIAB==0 & trainc$RF_GHYPE==0 &  trainc$RF_PPTERM==0, 0, 1)
   
   #-----------------------------------------------------------------------
   # Fit spline model.
   #-----------------------------------------------------------------------
   
   require(mgcv)
   
   # These are cubic p-splines, with 9 basis elements, with a second-order penalty.
   # See: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.terms.html
   fitmod = gam(cbind(LBW_yes, LBW_total) ~ s(combgest, bs='ps', m=c(2,2)) +
                  MRACE6 + mager + RF_GDIAB + RF_GHYPE + previs+ wtgain + CIG_REC + RF_PPTERM + SEX,
                data = trainc, family='binomial')
   
   # Predict on training set.
   predict_is = predict(fitmod, newdata=train, type='response', se=T)
   phat_is = as.vector(as.numeric(predict_is$fit)) 
   se_is = as.vector(as.numeric(predict_is$se.fit))
   
   # Predict on test set.
   predict_oos = predict(fitmod, newdata=test, type='response', se=T)
   phat_oos = as.vector(as.numeric(predict_oos$fit))
   se_oos = as.vector(as.numeric(predict_oos$se.fit))  
   
   fit = list('phat_is' = phat_is, 'se_is' = se_is, 'phat_oos' = phat_oos, 'se_oos' = se_oos)
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   return(list('fit' = fit,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos - 1.96 * se_oos,
               'phat_oos_ub' = phat_oos + 1.96 * se_oos,
               'binll_oos' = binll_oos))
}



########################################################################
# Random Forest fit utility.
########################################################################

rf_fit_util = function(train, test){
   
   require(ranger)
   
   #-------------------------------------------------------------------
   # Fit model to train data.set
   #-------------------------------------------------------------------
   rf_is = ranger(formula = factor(LBW) ~ combgest + MRACE6 + mager + RF_GDIAB + RF_GHYPE + previs+ wtgain + CIG_REC + RF_PPTERM + SEX, 
                  data=train, 
                     keep.inbag=T, oob.error=F, probability=T)
   
   # rf_is_pred = predict(rf_is, data=train, type='se')
   # phat_is = rf_is_pred$predictions
   # se_is = rf_is_pred$se
   
   #-----------------------------------------------------------------------
   # Predict on test set.
   #-----------------------------------------------------------------------
   test$LBW = factor(test$LBW)
   rf_oos_pred = predict(rf_is, data=test, type='se', probability=T)
   
   phat_oos = rf_oos_pred$predictions[,2]
   se_oos = rf_oos_pred$se[,2]
   
   #-------------------------------------------------------------------
   # Calculate fit info.
   #-------------------------------------------------------------------
   # Use test info only; don't include the patient panel obs.  This uses a global variable in dataset, testtrain.
   yobs_pred = as.numeric(test$LBW)
   binll_oos = binloglik(yobs_pred, phat_oos)
   
   return(list('fit' = rf_is,
               'phat_oos' = phat_oos,
               'phat_oos_lb' = phat_oos - 1.96 * se_oos,
               'phat_oos_ub' = phat_oos + 1.96 * se_oos,
               'phat_oos_se' = se_oos,
               'binll_oos' = binll_oos))
}




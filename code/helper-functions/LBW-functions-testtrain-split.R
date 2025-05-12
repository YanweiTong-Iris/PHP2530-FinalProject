########################################################################
# Split test/train data based on stratified sample
########################################################################

testtrain = function(data, test_pct = .3){
   
   ### Stratified sampling based on all combos of preterm/nd and diab_htn.
   
   # Library for left_join.
   require(dplyr)
   require(sampling)
   
   # Combine nd and preterm into one temp variable for convenience.
   #grps = expand.grid('preterm' = c(0,1), 
   #                   'diab_htn' = factor(c('Neither','Diabetes','Htn','Both'), 
   #                                       levels = c('Neither', 'Diabetes', 'Htn','Both')))
   
   grps = expand.grid('preterm' = c(0,1))
   
   grps$strata = 1:nrow(grps)
   
   # Match two columns to get strata group (out of possible 12 groups).
  # data$strata = suppressWarnings(left_join(data, grps, by=c('preterm','diab_htn'))$strata)
   data$strata = suppressWarnings(left_join(data, grps, by=c('preterm'))$strata)
   
   # Vectors to hold row indices of test and train datasets.
   test_rowidx = NULL
   train_rowidx = NULL
   
   # Loop through strata, and assign to train/test set.
   for(i in 1:length(grps$strata)){
      
      # Unique ids in strata i.
      ids_in_strata = unique(data$id[which(data$strata==i)])
      
      # Sample size of test set.
      sampsize = round(length(ids_in_strata) * test_pct)
      
      test_ids = sample(ids_in_strata, sampsize, replace=F)
      train_ids = ids_in_strata[-which(ids_in_strata %in% test_ids)]
      
      test_rowidx = c(test_rowidx, which(data$id %in% test_ids))
      train_rowidx = c(train_rowidx, which(data$id %in% train_ids))
      
   }
   
   # Test and train datasets.
   te = data[test_rowidx,]
   tr = data[train_rowidx,]
   
   # Return output
   return(list('test' = te, 'train' = tr))
}

########################################################################
# Upsampled data set for binll validation.
########################################################################

upsample = function(df, n=1000, sb_pct = .5){

   # Library for left_join.
   require(dplyr)
   require(sampling)

   # Combine nd and preterm into one temp variable for convenience.
   grps = expand.grid('combgest' = 30:42)
   grps$strata = 1:nrow(grps)

   # Match two columns to get strata group (out of possible 12 groups).
   df$strata = suppressWarnings(left_join(df, grps, by=c('combgest'))$strata)

   keep_preterm0 = NULL
   keep_preterm1 = NULL

   sampsizes = round(as.vector(table(df$strata)) / nrow(df) * n,0)
   sampsizes_preterm0 = ceiling(sampsizes * (1-sb_pct))
   sampsizes_preterm1 = ceiling(sampsizes * (sb_pct))

   # Loop through strata, and assign to train/test set.
   for(i in 1:length(grps$strata)){
     
     ids_preterm0 <- df$id[which(df$strata == i & df$preterm == 0)]
     ids_preterm1 <- df$id[which(df$strata == i & df$preterm == 1)]
     
     if (length(ids_preterm0) > 0 && sampsizes_preterm0[i] > 0) {
       samp_preterm0 <- sample(ids_preterm0, size = sampsizes_preterm0[i], replace = TRUE)
       keep_preterm0 <- c(keep_preterm0, samp_preterm0)
     }
     
     if (length(ids_preterm1) > 0 && sampsizes_preterm1[i] > 0) {
       samp_preterm1 <- sample(ids_preterm1, size = sampsizes_preterm1[i], replace = TRUE)
       keep_preterm1 <- c(keep_preterm1, samp_preterm1)
     }
   }

   keep_ids = c(keep_preterm0, keep_preterm1)
   keep_rows = which(df$id %in% keep_ids)

   # Test and train datasets.
   df = dplyr::left_join(data.frame('id'=keep_ids, drop=F), df, by='id', keep=FALSE) %>% 
     dplyr::select(-strata)

   # Return output
   return(df)
}

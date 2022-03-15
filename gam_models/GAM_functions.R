library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)

#FIT GAM SMOOTH FUNCTION
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit.smooth <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE){
  
#Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
#GAM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to 0

#GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  
  #Calculate the magnitude and significance of the smooth term effect by comparing full and reduced models
  ##Compare a full model GAM (with the smooth term) to a nested, reduced model (with covariates only)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, covariates)) #no smooth term
  gam.nullmodel <- gam(nullmodel, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ##Full versus reduced model anova p-value
  anova.smooth.pvalue <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  ### effect direction
  mean.derivative <- mean(derv$derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    partialRsq <- partialRsq*-1}
  
#Derivative-based temporal characteristics
  #Age of developmental change onset
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    change.onset <- min(derv$data[derv$sig==T])} #find first age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ #if gam derivative is never significant
    change.onset <- NA} #assign NA
  
  #Age of maximal developmental change
  if(sum(derv$sig) > 0){ 
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value significant derivatives
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    window.peak.change <- derv$data[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is greatest in absolute magnitude
    peak.change <- mean(window.peak.change)} #identify the age of peak developmental change
  if(sum(derv$sig) == 0){ 
    peak.change <- NA}  
  
  #Age of decrease onset
  if(sum(derv$sig) > 0){ 
    decreasing.range <- derv$data[derv$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
    if(length(decreasing.range) > 0)
      decrease.onset <- min(decreasing.range) #find youngest age with a significant negative derivative
    if(length(decreasing.range) == 0)
      decrease.onset <- NA}
  if(sum(derv$sig) == 0){
    decrease.onset <- NA}  
    
  #Age of increase offset
  if(sum(derv$sig) > 0){ 
    increasing.range <- derv$data[derv$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
    if(length(increasing.range) > 0)
      increase.offset <- max(increasing.range) #find oldest age with a significant positive derivative
    if(length(increasing.range) == 0)
      increase.offset <- NA}
  if(sum(derv$sig) == 0){ 
    increase.offset <- NA}  
  
  #Age of maturation
  if(sum(derv$sig) > 0){ 
    change.offset <- max(derv$data[derv$sig==T])} #find last age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ 
    change.offset <- NA}  
  
  full.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq, anova.smooth.pvalue, change.onset, peak.change, decrease.onset, increase.offset, change.offset)
  stats.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq, anova.smooth.pvalue)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}

#PREDICT GAM SMOOTH FUNCTION
##Function to predict fitted values of a measure based on a fitted GAM smooth (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and a prediction df
gam.smooth.predict <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, increments){

#Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
#Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 

#Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
#Generate predictions based on the gam model and predication data frame
  predicted.smooth <- fitted_values(object = gam.model, data = pred)
  predicted.smooth <- predicted.smooth %>% select(all_of(smooth_var), fitted, se, lower, upper)
  
#Determine the smooth_var index at which y is maximal, based on the predicted smooth
  maxy <- max(predicted.smooth$fitted)
  peak <- predicted.smooth[,smooth_var][predicted.smooth$fitted == maxy]
  
  smooth.fit <- list(parcel, peak, predicted.smooth)
  return(smooth.fit)
}

#POSTERIOR DISTRIBUTION SMOOTHS FUNCTION
##Function to simulate the posterior distribution from a fitted GAM, calculate smooths for individual posterior draws, and return smooth max and min values + 95% credible intervals
gam.posterior.smooths <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_draws = TRUE){

#Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior smooths estimated
  np <- as.numeric(increments) #number of smooth_var increments to predict fit at from posterior model coefficients
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  set.seed(10)

#Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
#Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
#Create a prediction data frame, used to estimate posterior model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to fit for each posterior draw
  for (v in c(1:length(theseVars))) { #fill the prediction df with data 
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
#Estimate posterior smooth functions (fitted values) from simulated GAM posterior distribution  
##Each of the posterior draws has a fitted spline (+ intercept + covariate coefficients) that includes the uncertainty in the estimated model coefficients  
##Akin to fitted_samples from gratia  
  Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
  sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance. 
  X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors that maps model parameters to the smooth fit (outcome measure scale)
  predicted.smooth.values <- X0 %*% t(sims) #generate posterior smooths (fitted y for each set of posterior draw model parameters)
  colnames(predicted.smooth.values) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
  predicted.smooth.values <- cbind(as.numeric(pred[,smooth_var]), predicted.smooth.values) #add smooth_var increments from pred df to first column
  colnames(predicted.smooth.values)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
  predicted.smooth.values <- as.data.frame(predicted.smooth.values) 

#Smooth minimum/maximum values and credible intervals
  #Smooth max + 95% credible interval
  max.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is largest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x){
                       round(predicted.smooth.values[,smooth_var][which.max(x)],2)
                     }))
  max.y.range <- t(max.y.range)
  max.y <- median(max.y.range) #median value 
  max.y.CI <- quantile(max.y.range, probs = c(0.025, 0.975)) #credible interval
  max.y.CI.lower <- max.y.CI[[1]]
  max.y.CI.upper <- max.y.CI[[2]]
  #Smooth min + 95% credible interval
  min.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is lowest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x){
                       round(predicted.smooth.values[,smooth_var][which.min(x)],2)
                     }))  
  min.y.range <- t(min.y.range)
  min.y <- median(min.y.range) #median value
  min.y.CI <- quantile(min.y.range, probs = c(0.025, 0.975)) #credible interval
  min.y.CI.lower <- min.y.CI[[1]]
  min.y.CI.upper <- min.y.CI[[2]]

  if(return_draws == TRUE)
    return(predicted.smooth.values)
  if(return_draws == FALSE)
    smooth.features <- list(parcel, max.y, max.y.CI.lower, max.y.CI.upper, min.y, min.y.CI.lower, min.y.CI.upper)
    names(smooth.features) <- c("parcel", sprintf("%s at max y", smooth_var), "max y credible interval lower", "max y credible interval upper", sprintf("%s at min y", smooth_var), "min y credible interval lower", "min y credible interval upper")
    return(smooth.features)
}

#DERIVATIVES FUNCTION
##Function to compute smooth derivatives for a main GAM model and for individual draws from the simulated posterior distribution
gam.derivatives <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_posterior_derivatives = TRUE){
  
#Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  set.seed(10)
  
#Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
#Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
#Create a prediction data frame, used to estimate posterior model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to fit for each posterior draw
  for (v in c(1:length(theseVars))) { #fill the prediction df with data 
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init) #prediction df
  pred2 <- pred #second prediction df
  pred2[,smooth_var] <- pred[,smooth_var] + EPS #finite differences
  
#Estimate smooth derivatives
  derivs <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = UNCONDITIONAL) #derivative at 200 indices of smooth_var with a simultaneous CI
  derivs.fulldf <- derivs %>% select(data, derivative, se, lower, upper)
  derivs.fulldf <- derivs.fulldf %>% mutate(significant = !(0 > lower & 0 < upper))
  derivs.fulldf$significant.derivative = derivs.fulldf$derivative*derivs.fulldf$significant
  colnames(derivs.fulldf) <- c(sprintf("%s", smooth_var), "derivative", "se", "lower", "upper", "significant", "significant.derivative")
  derivs <- derivs %>% select(data, derivative)
  colnames(derivs) <- c(sprintf("%s", smooth_var), "derivative")
  derivs[,smooth_var] <- round(derivs[,smooth_var], 3)
  derivs <- as.data.frame(derivs)
  
#Estimate posterior smooth derivatives from simulated GAM posterior distribution
  Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
  sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance. 
  X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors for pred
  X1 <- predict(gam.model, newdata = pred2, type = "lpmatrix") #get matrix of linear predictors for pred2
  Xp <- (X1 - X0) / EPS 
  posterior.derivs <- Xp %*% t(sims) #Xp * simulated model coefficients = simulated derivatives. Each column of posterior.derivs contains derivatives for a different draw from the simulated posterior distribution
  posterior.derivs <- as.data.frame(posterior.derivs)
  colnames(posterior.derivs) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
  posterior.derivs <- cbind(as.numeric(pred[,smooth_var]), posterior.derivs) #add smooth_var increments from pred df to first column
  colnames(posterior.derivs)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
  posterior.derivs[,smooth_var] <- round(posterior.derivs[,smooth_var], 3)
  
  if(return_posterior_derivatives == FALSE)
    return(derivs.fulldf)
  if(return_posterior_derivatives == TRUE)
    posterior.derivs <- merge(derivs, posterior.derivs, by="age")
    return(posterior.derivs)
}
  
#FIT GAM SMOOTH WITH A COVARIATE OF INTEREST FUNCTION
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariate of interest + control covariates)) and save out statistics for the first covariate
gam.fit.covariate <- function(measure, atlas, dataset, region, smooth_var, covariate.interest, covariates.noninterest, knots, set_fx = FALSE){

#Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s + %s", region, smooth_var, knots, set_fx, covariate.interest, covariates.noninterest))
  gam.model <- gam(modelformula, method = "REML", data=gam.data)
  gam.results <- summary(gam.model)
  
#GAM statistics
  #t-value for the covariate of interest term and GAM-based significance of this term
  gam.cov.tvalue <- gam.results$p.table[2,3]
  #GAM based significance of the term
  gam.cov.pvalue <- gam.results$p.table[2,4]
  
  #Calculate the magnitude and significance of the covariate of interest effect by comparing full and reduced models
  ##Compare a full model GAM (with the covariate of interst) to a nested, reduced model (without covariate of interst)
  nullmodel <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates.noninterest))
  gam.nullmodel <- gam(nullmodel, method = "REML", data=gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ##Full versus reduced model anova p-value
  anova.cov.pvalue <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.pvalue <- 1}
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  ### effect direction
  if(gam.cov.tvalue < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1}
  
  results <- cbind(parcel, gam.cov.tvalue, gam.cov.pvalue, anova.cov.pvalue, partialRsq)
  return(results)
}

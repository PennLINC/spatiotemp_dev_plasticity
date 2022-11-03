#Model Fitting: Developmental Environment-Dependent Variation in Regional Fluctuation Amplitude

library(dplyr)
library(cifti)
source("/cbica/projects/spatiotemp_dev_plasticity/code/spatiotemp_dev_plasticity/gam_models/GAM_functions.R")

# Read in Data
fluctuations.glasser.pnc <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/PNC/glasserfluctuations_demographics_finalsample.csv")
environment <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/sample_info/PNC/n1601_go1_environment_factor_scores_tymoore_20150909.csv")
environment <- environment %>% select(-bblid)
fluctuations.glasser.pnc <- merge(fluctuations.glasser.pnc, environment, by="scanid", sort = F)
fluctuations.glasser.pnc$sex <- as.factor(fluctuations.glasser.pnc$sex)

glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T) #glasser parcel names in order of surface data



#### Environment Factor Scores Across Age ####

cor.test(fluctuations.glasser.pnc$envSES, fluctuations.glasser.pnc$age) #t = 0.47644, df = 1031, p-value = 0.6339, cor = 0.01483646 



#### Fit GAM Models (Environment Effects) ####

## Controlling for age, sex, motion

gam.env.glasser <- matrix(data=NA, nrow=360, ncol=5) #matrix to save gam.fit.covariate output to

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$label[row] 
  GAM.RESULTS <- gam.fit.covariate(measure = "fluctuations", atlas = "glasser", dataset = "pnc", region = region, smooth_var = "age", covariate.interest = "envSES", covariates.noninterest = "sex + RMSmotion", knots = 3, set_fx = TRUE) #run the gam.fit.covariate function
  gam.env.glasser[row,] <- GAM.RESULTS} #and append results to output df 

gam.env.glasser <- as.data.frame(gam.env.glasser)
colnames(gam.env.glasser) <- c("label","GAM.env.tvalue","GAM.env.pvalue","Anova.env.pvalue","GAM.env.partialR2")
cols = c(2:5)    
gam.env.glasser[,cols] = apply(gam.env.glasser[,cols], 2, function(x) as.numeric(as.character(x)))

## Controlling for age, sex, motion, and parental education

fluctuations.glasser.pnc$parental.education <- fluctuations.glasser.pnc %>% select(medu1, fedu1) %>% rowMeans(na.rm = TRUE) #mean parental education measure
fluctuations.glasser.subset <- fluctuations.glasser.pnc %>% filter(parental.education != "NaN") #remove 7 participants without mother or father education data

gam.env.edu.glasser <- matrix(data=NA, nrow=360, ncol=5) #matrix to save gam.fit.covariate output to

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$label[row] 
  GAM.RESULTS <- gam.fit.covariate(measure = "fluctuations", atlas = "glasser", dataset = "subset", region = region, smooth_var = "age", covariate.interest = "envSES", covariates.noninterest = "sex + RMSmotion + parental.education", knots = 3, set_fx = TRUE) #run the gam.fit.covariate function
  gam.env.edu.glasser[row,] <- GAM.RESULTS} #and append results to output df 

gam.env.edu.glasser <- as.data.frame(gam.env.edu.glasser)
colnames(gam.env.edu.glasser) <- c("label","GAM.env.edu.tvalue","GAM.env.edu.pvalue","Anova.env.edu.pvalue","GAM.env.edu.partialR2")
cols = c(2:5)    
gam.env.edu.glasser[,cols] = apply(gam.env.edu.glasser[,cols], 2, function(x) as.numeric(as.character(x)))

## Save results
gam.env.glasser <- merge(gam.env.glasser, gam.env.edu.glasser, by="label", sort = F)
write.csv(gam.env.glasser, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/fluctuationamplitude_environment_statistics_glasser.csv", row.names = F, quote = F)



### Fit GAM Models (Environment Effects): Child, Adolescent, Early Adult Samples ####

## Child sample
fluctuations.glasser.child <- fluctuations.glasser.pnc %>% filter(age >= 8) %>% filter(age <= 12)
fluctuations.glasser.child$dev.stage <- c("child")

gam.env.child <- matrix(data=NA, nrow=360, ncol=5) #matrix to save gam.fit.covariate output to

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$label[row] 
  GAM.RESULTS <- gam.fit.covariate(measure = "fluctuations", atlas = "glasser", dataset = "child", region = region, smooth_var = "age", covariate.interest = "envSES", covariates.noninterest = "sex + RMSmotion", knots = 3, set_fx = TRUE) #run the gam.fit.covariate function
  gam.env.child[row,] <- GAM.RESULTS} #and append results to output df 

gam.env.child <- as.data.frame(gam.env.child)
colnames(gam.env.child) <- c("label","GAM.env.tvalue.child","GAM.env.pvalue.child","Anova.env.pvalue.child","GAM.env.partialR2.child")
cols = c(2:5)    
gam.env.child[,cols] = apply(gam.env.child[,cols], 2, function(x) as.numeric(as.character(x)))

## Adolescent sample
fluctuations.glasser.adolescent <- fluctuations.glasser.pnc %>% filter(age >= 13) %>% filter(age <= 17)
fluctuations.glasser.adolescent$dev.stage <- c("adolescent")

gam.env.adolescent <- matrix(data=NA, nrow=360, ncol=5) #matrix to save gam.fit.covariate output to

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$label[row] 
  GAM.RESULTS <- gam.fit.covariate(measure = "fluctuations", atlas = "glasser", dataset = "adolescent", region = region, smooth_var = "age", covariate.interest = "envSES", covariates.noninterest = "sex + RMSmotion", knots = 3, set_fx = TRUE) #run the gam.fit.covariate function
  gam.env.adolescent[row,] <- GAM.RESULTS} #and append results to output df 

gam.env.adolescent <- as.data.frame(gam.env.adolescent)
colnames(gam.env.adolescent) <- c("label","GAM.env.tvalue.adolescent","GAM.env.pvalue.adolescent","Anova.env.pvalue.adolescent","GAM.env.partialR2.adolescent")
cols = c(2:5)    
gam.env.adolescent[,cols] = apply(gam.env.adolescent[,cols], 2, function(x) as.numeric(as.character(x)))

## Early Adult sample
fluctuations.glasser.adult <- fluctuations.glasser.pnc %>% filter(age >= 18) %>% filter(age <= 23)
fluctuations.glasser.adult$dev.stage <- c("adult")

gam.env.adult <- matrix(data=NA, nrow=360, ncol=5) #matrix to save gam.fit.covariate output to

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$label[row] 
  GAM.RESULTS <- gam.fit.covariate(measure = "fluctuations", atlas = "glasser", dataset = "adult", region = region, smooth_var = "age", covariate.interest = "envSES", covariates.noninterest = "sex + RMSmotion", knots = 3, set_fx = TRUE) #run the gam.fit.covariate function
  gam.env.adult[row,] <- GAM.RESULTS} #and append results to output df 

gam.env.adult <- as.data.frame(gam.env.adult)
colnames(gam.env.adult) <- c("label","GAM.env.tvalue.adult","GAM.env.pvalue.adult","Anova.env.pvalue.adult","GAM.env.partialR2.adult")
cols = c(2:5)    
gam.env.adult[,cols] = apply(gam.env.adult[,cols], 2, function(x) as.numeric(as.character(x)))

## Save results
df.list <- list(gam.env.child, gam.env.adolescent, gam.env.adult)
gam.env.developmenetalstage <- Reduce(function(x,y) merge(x,y, all=TRUE, sort=F), df.list) 
write.csv(gam.env.developmenetalstage, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/fluctuationamplitude_environment_devstage_statistics_glasser.csv", row.names = F, quote = F)



#### Developmental Trajectories (Fitted Values) by Environment Factor Score Percentiles ####
np <- 200 #number of ages to predict fluctuation amplitude at

## 10th percentile neighborhood factor score
gam.smooths.lowenvSES <- matrix(data=NA, ncol=8) 
colnames(gam.smooths.lowenvSES) <- c("age","fitted","se","lower","upper","fitted.centered","index","label")

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$label[row] 
  SMOOTH.ESTIMATES <- gam.smooth.predict.covariateinteraction(measure = "fluctuations", atlas = "glasser", dataset = "pnc",
                                                              region = region, smooth_var = "age", int_var = "envSES", 
                                                              int_var.predict = quantile(fluctuations.glasser.pnc$envSES, c(.1)), 
                                                              covariates = "sex + RMSmotion", knots = 3, set_fx = TRUE, increments = np) #run gam.smooth.predict.covariateinteraction function
  SMOOTH.ESTIMATES$index <- rep(x=row, np) #region index
  SMOOTH.ESTIMATES$label <- rep(x=region, np) #label
  gam.smooths.lowenvSES <- rbind(gam.smooths.lowenvSES, SMOOTH.ESTIMATES)
}
gam.smooths.lowenvSES <- gam.smooths.lowenvSES[-1,] #remove empty initialization row

## 90th percentile neighborhood factor score
gam.smooths.highenvSES <- matrix(data=NA, ncol=8) 
colnames(gam.smooths.highenvSES) <- c("age","fitted","se","lower","upper","fitted.centered","index","label")

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$label[row] 
  SMOOTH.ESTIMATES <- gam.smooth.predict.covariateinteraction(measure = "fluctuations", atlas = "glasser", dataset = "pnc",
                                                              region = region, smooth_var = "age", int_var = "envSES", 
                                                              int_var.predict = quantile(fluctuations.glasser.pnc$envSES, c(.9)), 
                                                              covariates = "sex + RMSmotion", knots = 3, set_fx = TRUE, increments = np) #run gam.smooth.predict.covariateinteraction function
  SMOOTH.ESTIMATES$index <- rep(x=row, np) #region index
  SMOOTH.ESTIMATES$label <- rep(x=region, np) #label
  gam.smooths.highenvSES <- rbind(gam.smooths.highenvSES, SMOOTH.ESTIMATES)
}
gam.smooths.highenvSES <- gam.smooths.highenvSES[-1,] #remove empty initialization row

## Save results
gam.smooths.lowenvSES$envSES <- c("low")
gam.smooths.highenvSES$envSES <- c("high")
gam.smooths.envSES <- rbind(gam.smooths.lowenvSES, gam.smooths.highenvSES)
write.csv(gam.smooths.envSES,"/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/fluctuationamplitude_age_predictedfits_byenv.csv", row.names = F, quote = F) 



#### Age-specific Environment Effect Stratification along the Sensorimotor-Association Axis ####
np <- 200 #number of ages to get the derivative at
npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)

S.A.axis.cifti <- read_cifti("/cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSLRVertex/SensorimotorAssociation_Axis_parcellated/SensorimotorAssociation.Axis.Glasser360.pscalar.nii") #S-A_ArchetypalAxis repo
S.A.axis <- as.data.frame(cbind(rank(S.A.axis.cifti$data), names(S.A.axis.cifti$Parcel)))
colnames(S.A.axis) <- c("SA.axis","orig_parcelname")
S.A.axis <- merge(S.A.axis, glasser.parcel.labels, by="orig_parcelname", sort = F)
rm(S.A.axis.cifti)

SNR.mask <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_glasser360.csv")
S.A.axis <- merge(S.A.axis, SNR.mask, by = "label", sort=F)

compute_axis_correlation <- function(){ 
  # Get posterior effect estimates
  gam.enveffects.byage <- map_dfr(glasser.parcel.labels$label, 
                                     function(x){gam.varyingcoefficients(measure = "fluctuations", atlas = "glasser", dataset = "pnc",
                                                                         region = as.character(x), smooth_var = "age", int_var = "envSES", 
                                                                         covariates = "sex + RMSmotion", knots = 3, set_fx = TRUE, increments = np, 
                                                                         draws = 1000, return_posterior_coefficients = TRUE)}) #run gam.varyingcoefficients
  # Compute and save correlations with S-A axis
  gam.enveffects.byage <- left_join(S.A.axis, gam.enveffects.byage, by="label", sort=F) #assign axis rank to each label
  gam.enveffects.byage <- gam.enveffects.byage %>% filter(SNR.mask != 0) #include only high SNR parcels (N=336) in analyses
  corr_values <- gam.enveffects.byage %>%
    group_by(draw,age) %>%
    do(SAcorrelation = cor(as.numeric(.$SA.axis), as.numeric(.$envSES.slope), method=c("spearman"))) %>% 
    unnest(cols = c(SAcorrelation))
  corr_values.wide <- corr_values %>% pivot_wider(names_from = "draw", values_from = "SAcorrelation", names_sort = FALSE)
  corr_values.wide <- corr_values.wide %>% select(contains("draw"))
  rm(gam.enveffects.byage)
  gc()
  return(corr_values.wide)
}

enveffect.SAaxis.posteriorcor <- rerun(10, compute_axis_correlation()) %>% bind_cols()
write.table(enveffect.SAaxis.posteriorcor, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/SAaxis_environmenteffects_correlation_byage_glasser.csv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")

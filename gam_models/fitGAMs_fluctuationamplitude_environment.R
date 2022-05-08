#Model Fitting: Developmental Environment-Dependent Variation in Regional Fluctuation Amplitude

library(dplyr)
source("/cbica/projects/spatiotemp_dev_plasticity/code/spatiotemp_dev_plasticity/gam_models/GAM_functions.R")

# Read in Data
fluctuations.glasser.pnc <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/PNC/glasserfluctuations_demographics_finalsample.csv")
environment <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/sample_info/PNC/n1601_go1_environment_factor_scores_tymoore_20150909.csv")
environment <- environment %>% select(-bblid)
fluctuations.glasser.pnc <- merge(fluctuations.glasser.pnc, environment, by="scanid", sort = F)

glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T) #glasser parcel names in order of surface data

# Fit GAM Models (Environment Effects)

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




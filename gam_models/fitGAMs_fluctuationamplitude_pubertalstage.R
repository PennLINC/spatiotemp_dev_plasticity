# MODEL FITTING: EFFECTS OF PUBERTAL STAGE

library(dplyr)
source("/cbica/projects/spatiotemp_dev_plasticity/code/spatiotemp_dev_plasticity/gam_models/GAM_functions.R")

# Read in Data
glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T) #glasser parcel names in order of surface data

fluctuations.glasser.pnc <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/PNC/glasserfluctuations_demographics_finalsample.csv")
pubertal.data <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/sample_info/PNC/Tanner_menarche_age.csv")
pubertal.data <- pubertal.data %>% select(bblid, tanner_boy_6, tanner_girl_7)
fluctuations.glasser.puberty <- merge(fluctuations.glasser.pnc, pubertal.data, by="bblid", sort=F)
fluctuations.glasser.puberty$sex <- as.factor(fluctuations.glasser.puberty$sex)

# Tanner stage and pubertal stage variables
fluctuations.glasser.puberty$tanner_stage <- ifelse(fluctuations.glasser.puberty$sex == 1, fluctuations.glasser.puberty$tanner_boy_6, fluctuations.glasser.puberty$tanner_girl_7)
fluctuations.glasser.puberty <- fluctuations.glasser.puberty %>% filter(tanner_stage != "NA") #949 participants with self-reported pubertal status

##176 pre-pubertal, 283 mid-pubertal, 490 post-pubertal
fluctuations.glasser.puberty$pubertal_stage[fluctuations.glasser.puberty$tanner_stage == 1] <- 1
fluctuations.glasser.puberty$pubertal_stage[fluctuations.glasser.puberty$tanner_stage == 2] <- 1
fluctuations.glasser.puberty$pubertal_stage[fluctuations.glasser.puberty$tanner_stage == 3] <- 1
fluctuations.glasser.puberty$pubertal_stage[fluctuations.glasser.puberty$tanner_stage == 4] <- 2
fluctuations.glasser.puberty$pubertal_stage[fluctuations.glasser.puberty$tanner_stage == 5] <- 3

fluctuations.glasser.puberty$pubertal_stage.ordered <- ordered(fluctuations.glasser.puberty$pubertal_stage, levels = c(1,2,3)) #ordinal categorical variable

# Fit GAM Models (Pubertal Stage Effects)
gam.fit.puberty <- function(measure, atlas, dataset, region, smooth_var, covariate.interest, covariates.noninterest, knots, set_fx = FALSE){
  
  ##Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s + %s", region, smooth_var, knots, set_fx, covariate.interest, covariates.noninterest))
  gam.model <- gam(modelformula, method = "REML", data=gam.data)
  gam.results <- summary(gam.model)
  
  ##Linear and quadractic pubertal stage effects
  gam.cov.linear.tvalue <- gam.results$p.table[2,3] 
  gam.cov.quadratic.tvalue <- gam.results$p.table[3,3] 
  
  ##Anova significance of pubertal stage
  nullmodel <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates.noninterest))
  gam.nullmodel <- gam(nullmodel, method = "REML", data=gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  anova.cov.pvalue <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.pvalue <- 1}
  
  ##Full versus reduced model partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
 
  results <- cbind(parcel, gam.cov.linear.tvalue, gam.cov.quadratic.tvalue, partialRsq, anova.cov.pvalue)
  return(results)
}


gam.puberty <- matrix(data=NA, nrow=360, ncol=5)

for(row in c(1:nrow(glasser.parcel.labels))){ 
  region <- glasser.parcel.labels$label[row] 
  PUBERTY.RESULTS <- gam.fit.puberty(measure = "fluctuations", atlas = "glasser", dataset = "puberty", 
                                     region = region, smooth_var = "age", covariate.interest = "pubertal_stage.ordered", 
                                     covariates.noninterest = "sex + RMSmotion", knots = 3, set_fx = TRUE) 
  gam.puberty[row,] <- PUBERTY.RESULTS} 

gam.puberty <- as.data.frame(gam.puberty)
colnames(gam.puberty) <- c("label","GAM.puberty.linear.tvalue","GAM.puberty.quadratic.tvalue","GAM.puberty.partialR2", "Anova.puberty.pvalue")
cols = c(2:5)    
gam.puberty[,cols] = apply(gam.puberty[,cols], 2, function(x) as.numeric(as.character(x)))

write.csv(gam.puberty, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/fluctuationamplitude_puberty_statistics_glasser.csv", row.names = F, quote = F)


# Fit GAM Models (Age Effects controlling for Pubertal Stage)

gam.age.glasser <- matrix(data=NA, nrow=360, ncol=10) 

for(row in c(1:nrow(glasser.parcel.labels))){ 
  region <- glasser.parcel.labels$label[row] 
  AGE.RESULTS <- gam.fit.smooth(measure = "fluctuations", atlas = "glasser", dataset = "puberty", 
                                region = region, smooth_var = "age", covariates = "sex + RMSmotion + pubertal_stage.ordered", 
                                knots = 3, set_fx = TRUE, stats_only = FALSE) 
  gam.age.glasser[row,] <- AGE.RESULTS} 

gam.age.glasser <- as.data.frame(gam.age.glasser)
colnames(gam.age.glasser) <- c("label", #region name
                               "GAM.age.Fvalue", #GAM F-value for the age smooth term
                               "GAM.age.pvalue", #GAM p-value for the age smooth term
                               "GAM.age.partialR2", #partial Rsq from age and age-null models
                               "Anova.age.pvalue", #Anova p-value comparing age and age-null models
                               "age.onsetchange", #age at which fluctuation amplitude starts significantly changing (first significant derivative)
                               "age.peakchange", #age at which fluctuation amplitude exhibits maximal change (largest significant derivative)
                               "minage.decrease", #age at which fluctuation amplitude starts significantly decreasing (first significant negative derivative)
                               "maxage.increase", #age at which fluctuation amplitude stops significantly increasing (last significant positive derivative)
                               "age.maturation") #age at which fluctuation amplitude stops changing (last significant derivative)
cols = c(2:10)    
gam.age.glasser[,cols] = apply(gam.age.glasser[,cols], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(gam.age.glasser, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/fluctuationamplitude_age_pubertycontrolled_statistics_glasser.csv", row.names = F, quote = F)

  
  
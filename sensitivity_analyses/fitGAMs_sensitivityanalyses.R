library(dplyr)
library(cifti)
library(purrr)
source("/cbica/projects/spatiotemp_dev_plasticity/code/spatiotemp_dev_plasticity/gam_models/GAM_functions.R")
np <- 200
npd <- 100

glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T)
fluctuations.glasser.pnc <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/PNC/glasserfluctuations_demographics_finalsample.csv")
fluctuations.glasser.pnc$sex <- as.factor(fluctuations.glasser.pnc$sex)

S.A.axis.cifti <- read_cifti("/cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSLRVertex/SensorimotorAssociation_Axis_parcellated/SensorimotorAssociation.Axis.Glasser360.pscalar.nii") #S-A_ArchetypalAxis repo
S.A.axis <- as.data.frame(cbind(rank(S.A.axis.cifti$data), names(S.A.axis.cifti$Parcel)))
colnames(S.A.axis) <- c("SA.axis","orig_parcelname")
S.A.axis <- merge(S.A.axis, glasser.parcel.labels, by="orig_parcelname", sort = F)
rm(S.A.axis.cifti)

SNR.mask <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_glasser360.csv")
S.A.axis <- merge(S.A.axis, SNR.mask, by = "label", sort=F)

compute_axis_correlation <- function(sensitivity.measure, sensitivity.dataset, sensitivity.type){
  # Get posterior derivatives
  if(sensitivity.type == "motion" || sensitivity.type ==  "psych" || sensitivity.type == "meannorm"){
  gam.derivatives.glasser <- map_dfr(glasser.parcel.labels$label, 
                                     function(x){gam.derivatives(measure = sensitivity.measure, atlas = "glasser", dataset = sensitivity.dataset, 
                                                                 region = as.character(x), smooth_var = "age", covariates = "sex + RMSmotion", 
                                                                 knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_derivatives = TRUE)}) 
  }
  if(sensitivity.type == "vascular"){
    gam.derivatives.glasser <- map_dfr(glasser.parcel.labels$label, 
                                       function(x){gam.derivatives(measure = sensitivity.measure, atlas = "glasser", dataset = sensitivity.dataset, 
                                                                   region = as.character(x), smooth_var = "age", covariates = sprintf("sex + RMSmotion + CBF_%s", str_replace(as.character(x), "-", ".")), 
                                                                   knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_derivatives = TRUE)}) 
  }
  if(sensitivity.type == "BOLD"){
    gam.derivatives.glasser <- map_dfr(glasser.parcel.labels$label, 
                                       function(x){gam.derivatives(measure = sensitivity.measure, atlas = "glasser", dataset = sensitivity.dataset, 
                                                                   region = as.character(x), smooth_var = "age", covariates = sprintf("sex + RMSmotion + BOLD_%s", str_replace(as.character(x), "-", ".")), 
                                                                   knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_derivatives = TRUE)}) 
  }
  
  # Compute and save correlations with S-A axis
  gam.derivatives.glasser <- left_join(S.A.axis, gam.derivatives.glasser, by="label", sort=F) #assign axis rank to each label
  gam.derivatives.glasser <- gam.derivatives.glasser %>% filter(SNR.mask != 0) #include only high SNR parcels (N=336) in analyses
  corr_values <- gam.derivatives.glasser %>%
    group_by(draw,age) %>%
    do(SAcorrelation = cor(as.numeric(.$SA.axis), as.numeric(.$posterior.derivative), method=c("spearman"))) %>%
    unnest(cols = c(SAcorrelation))
  corr_values.wide <- corr_values %>% pivot_wider(names_from = "draw", values_from = "SAcorrelation", names_sort = FALSE)
  corr_values.wide <- corr_values.wide %>% select(contains("draw"))
  rm(gam.derivatives.glasser)
  gc()
  return(corr_values.wide)
}

#### LOW MOTION SAMPLE ####
#remove 343 individuals from the main sample, low motion sample N = 690

fluctuations.glasser.lowmotion <- fluctuations.glasser.pnc %>% filter(RMSmotion < 0.075) 

#stats
statistics.lowmotion <- map_dfr(glasser.parcel.labels$label, 
                        function(x){as.data.frame(gam.fit.smooth(measure = "fluctuations", atlas = "glasser", dataset = "lowmotion", 
                                   region = as.character(x), smooth_var = "age", covariates = "sex + RMSmotion", 
                                   knots = 3, set_fx = TRUE, stats_only = FALSE))}) %>% 
                        set_names(c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.partialR2","Anova.age.pvalue","age.onsetchange","age.peakchange","minage.decrease","maxage.increase","age.maturation"))
statistics.lowmotion[,2:10] = apply(statistics.lowmotion[,2:10], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(statistics.lowmotion, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/lowmotion/statistics_lowmotion.csv", row.names = F, quote = F)
rm(statistics.lowmotion)

#smooths
smooths.lowmotion <- map_dfr(glasser.parcel.labels$label, 
                     function(x){as.data.frame(gam.estimate.smooth(measure = "fluctuations", atlas = "glasser", dataset = "lowmotion", 
                                                      region = as.character(x), smooth_var = "age", covariates = "sex + RMSmotion", 
                                                      knots = 3, set_fx = TRUE, increments = np))}) 

regionlist <- map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(as.character(x), np))}) %>% set_names(c("region"))
indexlist <-  map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(which(glasser.parcel.labels$label == as.character(x), arr.ind=TRUE), np))}) %>% set_names(c("index"))      
smooths.lowmotion <- cbind(smooths.lowmotion, indexlist, regionlist)
write.csv(smooths.lowmotion, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/lowmotion/smooths_lowmotion.csv", row.names = F, quote = F)
rm(smooths.lowmotion)

#derivative correlation
axiscorrelation.lowmotion <- rerun(10, compute_axis_correlation(sensitivity.measure = "fluctuations", sensitivity.dataset = "lowmotion", sensitivity.type = "motion")) %>% bind_cols()
write.table(axiscorrelation.lowmotion, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/lowmotion/axiscorrelation_lowmotion.csv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
rm(axiscorrelation.lowmotion)
rm(fluctuations.glasser.lowmotion)
gc()


#### PSYCHIATRIC EXCLUDES (LTN SAMPLE) ####
#remove 140 individuals from the main sample with psych hospitalization or current med use, psychiatric exclusion N = 893

LTNexclude <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/sample_info/PNC/n1601_health_20170421.csv")
fluctuations.glasser.ltn <- left_join(fluctuations.glasser.pnc, LTNexclude, by = "scanid", sort = F)
fluctuations.glasser.ltn <- fluctuations.glasser.ltn %>% filter(ltnExcludev2 == 0) 

#stats
statistics.ltn <- map_dfr(glasser.parcel.labels$label, 
                  function(x){as.data.frame(gam.fit.smooth(measure = "fluctuations", atlas = "glasser", dataset = "ltn", 
                                                                         region = as.character(x), smooth_var = "age", covariates = "sex + RMSmotion", 
                                                                         knots = 3, set_fx = TRUE, stats_only = FALSE))}) %>% 
                  set_names(c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.partialR2","Anova.age.pvalue","age.onsetchange","age.peakchange","minage.decrease","maxage.increase","age.maturation"))
statistics.ltn[,2:10] = apply(statistics.ltn[,2:10], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(statistics.ltn, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/ltn/statistics_ltn.csv", row.names = F, quote = F)
rm(statistics.ltn)

#smooths
smooths.ltn <- map_dfr(glasser.parcel.labels$label, 
                             function(x){as.data.frame(gam.estimate.smooth(measure = "fluctuations", atlas = "glasser", dataset = "ltn", 
                                                                           region = as.character(x), smooth_var = "age", covariates = "sex + RMSmotion", 
                                                                           knots = 3, set_fx = TRUE, increments = np))}) 

regionlist <- map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(as.character(x), np))}) %>% set_names(c("region"))
indexlist <-  map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(which(glasser.parcel.labels$label == as.character(x), arr.ind=TRUE), np))}) %>% set_names(c("index"))      
smooths.ltn <- cbind(smooths.ltn, indexlist, regionlist)
write.csv(smooths.ltn, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/ltn/smooths_ltn.csv", row.names = F, quote = F)
rm(smooths.ltn)

#derivative correlation
axiscorrelation.ltn <- rerun(10, compute_axis_correlation(sensitivity.measure = "fluctuations", sensitivity.dataset = "ltn", sensitivity.type = "psych")) %>% bind_cols()
write.table(axiscorrelation.ltn, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/ltn/axiscorrelation_ltn.csv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
rm(axiscorrelation.ltn)
rm(fluctuations.glasser.ltn)
gc()


#### CBF Control ####
#control regional fluctuation amplitude GAMs for regional cerebral blood flow (main sample with CBF data, N = 1002)

cbf.glasser.pnc <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/CBF/CBF_PNC_Glasser360.csv")

regionnames <- paste0("CBF_",colnames(cbf.glasser.pnc[1:360]))
demoname <- c("rbcid")
names <- c(regionnames, demoname)
colnames(cbf.glasser.pnc) <- names
fluctuations.glasser.cbf <- merge(fluctuations.glasser.pnc, cbf.glasser.pnc, by="rbcid", sort = F)

#stats
statistics.cbf <- map_dfr(glasser.parcel.labels$label, 
                          function(x){as.data.frame(gam.fit.smooth(measure = "fluctuations", atlas = "glasser", dataset = "cbf", 
                                                                   region = as.character(x), smooth_var = "age", covariates = sprintf("sex + RMSmotion + CBF_%s", str_replace(as.character(x), "-", ".")), 
                                                                   knots = 3, set_fx = TRUE, stats_only = FALSE))}) %>% 
  set_names(c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.partialR2","Anova.age.pvalue","age.onsetchange","age.peakchange","minage.decrease","maxage.increase","age.maturation"))
statistics.cbf[,2:10] = apply(statistics.cbf[,2:10], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(statistics.cbf, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/cbf/statistics_cbf.csv", row.names = F, quote = F)
rm(statistics.cbf)

#smooths
smooths.cbf <- map_dfr(glasser.parcel.labels$label, 
                       function(x){as.data.frame(gam.estimate.smooth(measure = "fluctuations", atlas = "glasser", dataset = "cbf", 
                                                                     region = as.character(x), smooth_var = "age", covariates = sprintf("sex + RMSmotion + CBF_%s", str_replace(as.character(x), "-", ".")), 
                                                                     knots = 3, set_fx = TRUE, increments = np))}) 

regionlist <- map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(as.character(x), np))}) %>% set_names(c("region"))
indexlist <-  map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(which(glasser.parcel.labels$label == as.character(x), arr.ind=TRUE), np))}) %>% set_names(c("index"))      
smooths.cbf <- cbind(smooths.cbf, indexlist, regionlist)
write.csv(smooths.cbf, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/cbf/smooths_cbf.csv", row.names = F, quote = F)
rm(smooths.cbf)

#derivative correlation
axiscorrelation.cbf <- rerun(10, compute_axis_correlation(sensitivity.measure = "fluctuations", sensitivity.dataset = "cbf", sensitivity.type = "vascular")) %>% bind_cols()
write.table(axiscorrelation.cbf, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/cbf/axiscorrelation_cbf.csv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
rm(axiscorrelation.cbf)
rm(fluctuations.glasser.cbf)
rm(cbf.glasser.pnc)
gc()


#### T2* Control ####
#control regional fluctuation amplitude GAMs for regional mean T2* signal (main sample, N = 1033)
#parcellated mean T2* maps were generated for each participant with /sensitivity_analyses/parcellated_meanT2star.R in study repo

bold.glasser.pnc <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Timeseries/fmriprep/meanBOLD_subxparcel_glasser360.csv")

regionnames <- paste0("BOLD_",colnames(bold.glasser.pnc[2:361]))
demoname <- c("rbcid")
names <- c(demoname,regionnames)
colnames(bold.glasser.pnc) <- names
fluctuations.glasser.T2 <- merge(fluctuations.glasser.pnc, bold.glasser.pnc, by="rbcid", sort = F)

#stats
statistics.T2 <- map_dfr(glasser.parcel.labels$label, 
                          function(x){as.data.frame(gam.fit.smooth(measure = "fluctuations", atlas = "glasser", dataset = "T2", 
                                                                   region = as.character(x), smooth_var = "age", covariates = sprintf("sex + RMSmotion + BOLD_%s", str_replace(as.character(x), "-", ".")), 
                                                                   knots = 3, set_fx = TRUE, stats_only = FALSE))}) %>% 
  set_names(c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.partialR2","Anova.age.pvalue","age.onsetchange","age.peakchange","minage.decrease","maxage.increase","age.maturation"))
statistics.T2[,2:10] = apply(statistics.T2[,2:10], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(statistics.T2, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/T2star/statistics_T2.csv", row.names = F, quote = F)
rm(statistics.T2)

#smooths
smooths.T2 <- map_dfr(glasser.parcel.labels$label, 
                       function(x){as.data.frame(gam.estimate.smooth(measure = "fluctuations", atlas = "glasser", dataset = "T2", 
                                                                     region = as.character(x), smooth_var = "age", covariates = sprintf("sex + RMSmotion + BOLD_%s", str_replace(as.character(x), "-", ".")), 
                                                                     knots = 3, set_fx = TRUE, increments = np))}) 

regionlist <- map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(as.character(x), np))}) %>% set_names(c("region"))
indexlist <-  map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(which(glasser.parcel.labels$label == as.character(x), arr.ind=TRUE), np))}) %>% set_names(c("index"))      
smooths.T2 <- cbind(smooths.T2, indexlist, regionlist)
write.csv(smooths.T2, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/T2star/smooths_T2.csv", row.names = F, quote = F)
rm(smooths.T2)

#derivative correlation
axiscorrelation.T2 <- rerun(10, compute_axis_correlation(sensitivity.measure = "fluctuations", sensitivity.dataset = "T2", sensitivity.type = "BOLD")) %>% bind_cols()
write.table(axiscorrelation.T2, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/T2star/axiscorrelation_T2.csv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
rm(axiscorrelation.T2)
gc()


#### Mean Normalized Fluctuation Amplitude ####
#divide fluctuation amplitude in each region by whole brain mean fluctuation amplitude

fluctuations.glasser.norm <- fluctuations.glasser.pnc[,2:361]
fluctuations.glasser.norm <- fluctuations.glasser.norm/rowMeans(fluctuations.glasser.norm) #parcel amplitude divided by global mean fluctuation amplitude)
fluctuations.glasser.norm$rbcid <- fluctuations.glasser.pnc$rbcid
fluctuations.glasser.norm$age <- fluctuations.glasser.pnc$age
fluctuations.glasser.norm$sex <- fluctuations.glasser.pnc$sex
fluctuations.glasser.norm$RMSmotion <- fluctuations.glasser.pnc$RMSmotion

#stats
statistics.norm <- map_dfr(glasser.parcel.labels$label, 
                          function(x){as.data.frame(gam.fit.smooth(measure = "fluctuations", atlas = "glasser", dataset = "norm", 
                                                                   region = as.character(x), smooth_var = "age", covariates = "sex + RMSmotion", 
                                                                   knots = 3, set_fx = TRUE, stats_only = FALSE))}) %>% 
  set_names(c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.partialR2","Anova.age.pvalue","age.onsetchange","age.peakchange","minage.decrease","maxage.increase","age.maturation"))
statistics.norm[,2:10] = apply(statistics.norm[,2:10], 2, function(x) as.numeric(as.character(x))) #format as numeric
write.csv(statistics.norm, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/norm/statistics_norm.csv", row.names = F, quote = F)
rm(statistics.norm)

#smooths
smooths.norm <- map_dfr(glasser.parcel.labels$label, 
                       function(x){as.data.frame(gam.estimate.smooth(measure = "fluctuations", atlas = "glasser", dataset = "norm", 
                                                                     region = as.character(x), smooth_var = "age", covariates = "sex + RMSmotion", 
                                                                     knots = 3, set_fx = TRUE, increments = np))}) 

regionlist <- map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(as.character(x), np))}) %>% set_names(c("region"))
indexlist <-  map_dfr(glasser.parcel.labels$label, function(x){as.data.frame(rep(which(glasser.parcel.labels$label == as.character(x), arr.ind=TRUE), np))}) %>% set_names(c("index"))      
smooths.norm <- cbind(smooths.norm, indexlist, regionlist)
write.csv(smooths.norm, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/norm/smooths_norm.csv", row.names = F, quote = F)
rm(smooths.norm)

#derivative correlation
axiscorrelation.norm <- rerun(10, compute_axis_correlation(sensitivity.measure = "fluctuations", sensitivity.dataset = "norm", sensitivity.type = "meannorm")) %>% bind_cols()
write.table(axiscorrelation.norm, "/cbica/projects/spatiotemp_dev_plasticity/FluctuationAmplitude/GAMRESULTS/sensitivity_analyses/norm/axiscorrelation_norm.csv", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
rm(axiscorrelation.norm)
rm(fluctuations.glasser.norm)
gc()
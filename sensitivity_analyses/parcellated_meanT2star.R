library(cifti)

participants <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/sample_info/PNC/PNC_timeseries_processed.txt", header=F) 
colnames(participants) <- c("rbcid") 

glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T)

# Compute parcelwise resting state functional MRI scan average T2*

compute_meanbold <- function(rbcid, atlas){
  if(atlas == "glasser360"){
    if(!file.exists(sprintf("/cbica/projects/spatiotemp_dev_plasticity/Timeseries/fmriprep/%1$s_fsLR_meanBOLD_glasser360.pscalar.csv", rbcid))){
      parcellated.bold.cifti <- read_cifti(sprintf("/cbica/projects/spatiotemp_dev_plasticity/Timeseries/fmriprep/%1$s_ses-PNC1_task-rest_acq-singleband_space-fsLR_atlas-Glasser_fmriprep_bold.ptseries.nii", rbcid))
      parcellated.meanbold <- rowMeans(parcellated.bold.cifti$data)
      parcellated.meanbold <- as.data.frame(parcellated.meanbold)
      parcellated.meanbold$parcelname <- names(parcellated.bold.cifti$Parcel)
      parcellated.meanbold$label <- glasser.parcel.labels$label
      write.csv(parcellated.meanbold, sprintf("/cbica/projects/spatiotemp_dev_plasticity/Timeseries/fmriprep/%1$s_fsLR_meanBOLD_glasser360.pscalar.csv", rbcid), quote=F, row.names=F)}}
}

for(sub in c(1:nrow(participants))){
  rbcid <- participants[sub,1]
  compute_meanbold(rbcid, "glasser360")
}

# Create Participant x Region mean BOLD spreadsheets

bold.subxparcel.matrix.glasser <- matrix(data = NA, nrow = nrow(participants), ncol = 361)
regionheaders <- glasser.parcel.labels$label
demoheaders <- c("rbcid")
colheaders <- as.matrix(c(demoheaders,regionheaders))
colnames(bold.subxparcel.matrix.glasser) <- colheaders

for(sub in c(1:nrow(participants))){
  rbcid=participants[sub,1]
  meanbold.data <- read.csv(sprintf("/cbica/projects/spatiotemp_dev_plasticity/Timeseries/fmriprep/%1$s_fsLR_meanBOLD_glasser360.pscalar.csv", rbcid))
  meanbold.data <- meanbold.data[,1]
  bold.subxparcel.matrix.glasser[sub,] <- cbind(rbcid, t(meanbold.data))
}

write.csv(bold.subxparcel.matrix.glasser, file = "/cbica/projects/spatiotemp_dev_plasticity/Timeseries/fmriprep/meanBOLD_subxparcel_glasser360.csv", row.names = F, quote = F)

library(cifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Users/valeriesydnor/Software/workbench/')
library(gifti)

#Resample vertex-wise SNR map from fsaverage5 surface mesh to fslr surface mesh
command1 = "-metric-resample /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fsaverage5_lefthemisphere.func.gii /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSaverage5/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSLRVertex/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fslr_lefthemisphere.func.gii -area-metrics /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSaverage5/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSLRVertex/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii"
ciftiTools::run_wb_cmd(command1, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)
command2 = "-metric-resample /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fsaverage5_righthemisphere.func.gii /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSaverage5/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSLRVertex/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii ADAP_BARY_AREA /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fslr_righthemisphere.func.gii -area-metrics /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSaverage5/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii /cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/FSLRVertex/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii"
ciftiTools::run_wb_cmd(command2, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

#Read in, binarize, and combine left and right vertex-wise SNR maps 
SNR.mask.lh <- read_gifti("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fslr_lefthemisphere.func.gii")
SNR.mask.lh$data$normal[SNR.mask.lh$data$normal > 0] <- 1 #binarize into mask, with low SNR vertices == 0
SNR.mask.rh <- read_gifti("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fslr_righthemisphere.func.gii")
SNR.mask.rh$data$normal[SNR.mask.rh$data$normal > 0] <- 1 #binarize into mask, with low SNR vertices == 0
SNR.mask <- as_cifti(cortexL = SNR.mask.lh$data$normal, cortexR = SNR.mask.rh$data$normal) #merge left and right vertex-wise SNR masks into one cifti
surfL_fname <- read_surf("/cbica/projects/spatiotemp_dev_plasticity/software/workbench/workbench_files/Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii")
surfR_fname <- read_surf("/cbica/projects/spatiotemp_dev_plasticity/software/workbench/workbench_files/Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii")
SNR.mask <- add_surf(SNR.mask, surfL=surfL_fname, surfR=surfR_fname)
#view_cifti_surface(SNR.mask) #visualize vertex-wise SNR mask
write_cifti(SNR.mask, "/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fslr.dscalar.nii")

#Parcellate vertex-wise SNR mask with study atlases 
command1 = "-cifti-parcellate /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fslr.dscalar.nii /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser_space-fsLR_den-32k_desc-atlas.dlabel.nii COLUMN /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_glasser360.pscalar.nii"
ciftiTools::run_wb_cmd(command1, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)
command2 = "-cifti-parcellate /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_fslr.dscalar.nii /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/Schaefer2018_400Parcels_17Networks_order.dlabel.nii COLUMN /cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_schaefer400.pscalar.nii"
ciftiTools::run_wb_cmd(command2, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

#Read in parcellated SNR masks, exclude parcels where >= 25% of vertices are low SNR (0), and save parcellated SNR mask csvs
SNR.mask.glasser <- cifti::read_cifti("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_glasser360.pscalar.nii")
SNR.mask.glasser <- as.data.frame(SNR.mask.glasser$data)
colnames(SNR.mask.glasser) <- c("SNR.mask")
glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T)
SNR.mask.glasser$label <- glasser.parcel.labels$label
SNR.mask.glasser$SNR.mask[SNR.mask.glasser$SNR.mask <= 0.75] <- 0 #exclude parcels where >= 25% of vertices have low SNR
SNR.mask.glasser$SNR.mask[SNR.mask.glasser$SNR.mask > 0.75] <- 1 #retain parcels where >75% of vertices have high SNR
sum(SNR.mask.glasser$SNR.mask == 0) #24 glasser parcels will be removed from analysis
write.csv(SNR.mask.glasser, "/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_glasser360.csv", quote = F, row.names = F)

SNR.mask.schaefer <- cifti::read_cifti("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_schaefer400.pscalar.nii")
SNR.mask.schaefer <- as.data.frame(SNR.mask.schaefer$data)
colnames(SNR.mask.schaefer) <- c("SNR.mask")
schaefer.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/schaefer400_regionlist.csv", header = T)
SNR.mask.schaefer$label <- schaefer.parcel.labels$label
SNR.mask.schaefer$SNR.mask[SNR.mask.schaefer$SNR.mask <= 0.75] <- 0
SNR.mask.schaefer$SNR.mask[SNR.mask.schaefer$SNR.mask > 0.75] <- 1
sum(SNR.mask.schaefer$SNR.mask == 0) #24 schaefer parcels will be removed from analysis
write.csv(SNR.mask.schaefer, "/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_schaefer400.csv", quote = F, row.names = F)

<br>
<br>
# Spontaneous Activity Development Unfolds Hierarchically Along the Sensorimotor-Association Axis 
*Spontaneous activity critically refines cortical circuits in the developing brain. In the rodent brain, spontaneous activity in sensory regions evolves from strong and synchronized to sparse and decorrelated as the cortex transitions from plastic to mature. Synchronized bursts in intrinsic activity are therefore considered a functional hallmark of plastic cortices. Here, we leverage this functional hallmark to characterize the developmental unfolding of activity-indexed plasticity in the human brain during youth. We hypothesized that age-related changes in spontaneous activity would be spatially and temporally governed by the sensorimotor-association (S-A) axis of cortical organization and provide evidence for hierarchical neurodevelopment in childhood and adolescence. Ultimately, elucidating spatially-localized windows of enhanced and diminished cortical plasticity will help to identify windows wherein youth educational and psychiatric interventions may be maximally beneficial.* 

### Project Lead
Valerie J. Sydnor

### Faculty Lead
Theodore D. Satterthwaite

### Analytic Replicator
Bart Larsen

### Collaborators 
Bart Larsen, Azeez Adebimpe, Maxwell A. Bertolero, Matthew Cieslak, Sydney Covitz, Yong Fan, Raquel E. Gur, Ruben C. Gur, David R. Roalf, Russell T. Shinohara, Dani S. Bassett, Theodore D. Satterthwaite

### Project Start Date
June 2021

### Current Project Status
Manuscript in preparation

### Datasets
RBC PNC-Health Exclude (primary) and LTN (sensitivity)

### Github Repository
<https://github.com/PennLINC/spatiotemp_dev_plasticity>

### Conference Presentations
- Oral presentation *Developmental Refinement of Spontaneous Activity Varies Across Sensorimotor and Association Cortices* selected for "From the Womb to Wisdom" oral session at The Organization for Human Brain Mapping Annual Meeting, June 2022. 
- Poster accepted at The Society of Biological Psychiatry Annual Meeting, April 2022 as a Travel Awardee.

### Cubic Project Directory
/cbica/projects/spatiotemp_dev_plasticity

```
CBF/: parcel-wise cerebral blood flow maps for each participant, generated with ASLPrep  
code/: directory with the spatiotemp_dev_plasticity github repo  
FluctuationAmplitude/PNC/: vertex-wise and parcel-wise functional amplitude maps for each participant, generated with xcp-d and connectome workbench  
FluctuationAmplitude/GAMRESULTS/: gam model outputs (effect sizes, p-values, fitted values, smooth estimates, smooth characteristics, derivatives)  
Maps/: surface parcellation files and SNR masks (parcellations directory) and S-A axis github repo (S-A_ArchetypalAxis directory)
Myelin/: myelin development maps including the age effect size map (r2), the age of maximal growth map (age of max slope), and the annualized rate of change map (annualized roc) from Baum et al., 2021  
sample_info/: sample demographics, factor scores, rbcid-bblid key, and final project participant list (PNC_FinalSample_N1033.csv)
software/: software used  
Structural/: freesurfer output for each participant  
Timeseries/: vertex-wise and parcel-wise fully processed resting fMRI timeseries data for each participant, generated with fmriprep and xcp-d  
```



<br>
<br>
# CODE DOCUMENTATION
The entire analytic workflow implemented in this project is described in the following sections with links to the corresponding github code. The workflow includes quantification of regional fluctuation amplitude, PNC sample selection, fitting of generalized additive models (GAMs), and characterization of relationships between fluctuation amplitude, age, environmental variability, and the sensorimotor-association axis. Scripts were implemented in the order outlined below.
<br>
### Fluctuation Amplitude Quantification
Resting state functional MRI data were processed with [fmriprep 20.2.3](https://hub.docker.com/layers/fmriprep/nipreps/fmriprep/20.2.3/images/sha256-102db5fe8b0a34298f2eb2fd5962ad99ff0a948d258cbf08736fcc1b845cad9f?context=explore) and [xcp-d 0.0.4](https://hub.docker.com/layers/xcp_abcd/pennlinc/xcp_abcd/0.0.4/images/sha256-317160b8078cf7978eaf9db6fef32df78864232cb8a8759a354832813d1faf02?context=explore) to quantify fluctuation amplitude at each vertex on the fslr 32k cortical surface. 

fmriprep 20.2.3 was run with the following parameters:

```bash
$ singularity run pennlinc-containers/.datalad/environments/fmriprep-20-2-3/image inputs/data prep participant --output-spaces MNI152NLin6Asym:res-2 --participant-label "$subid" --force-bbr --cifti-output 91k -v -v
```

xcp-d  0.0.4 was run with the following parameters: 

```bash
$ singularity run pennlinc-containers/.datalad/environments/xcp-abcd-0-0-4/image inputs/data/fmriprep xcp participant --despike --lower-bpf 0.01 --upper-bpf 0.08 --participant_label $subid -p 36P -f 10 –cifti
```

Vertex-wise fluctuation amplitude maps were then parcellated with [/fluctuation_amplitude/parcellate.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/fluctuation_amplitude/parcellate.Rmd) to quantify mean fluctuation amplitude in each cortical region. Regions were defined with the [HCP multimodal atlas](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/glasser_space-fsLR_den-32k_desc-atlas.dlabel.nii) (i.e. Glasser 360, primary analyses) and with the [Schaefer 400 atlas](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/Schaefer2018_400Parcels_17Networks_order.dlabel.nii) (sensitivity analysis). 

Associations between regional fluctuation amplitude and demographics were only examined in brain regions exhibiting high signal to noise ratio (SNR) in PNC functional MRI data. A vertex level SNR map was generated as in Cui et al., 2020, Neuron and parcellated with Glasser 360 and Schaefer 400 atlases with the script [/fluctuation_amplitude/SNR_mask.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/fluctuation_amplitude/SNR_mask.R). Regions wherein >= 25% of vertices had attenuated signal were excluded from analyses.


### Sample Construction
The final study sample was constructed with [/sample_construction/finalsample.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/sample_construction/finalsample.Rmd). The final sample was generated from the 1374 PNC participants with dominant group ses-PNC1_task-rest_acq-singleband data (non-variant CuBIDS acquisitions). The following exclusions were applied to generate the final sample (N = 1033):

> Health exclude: 120 participants with medical problems that could impact brain function or incidentally-encountered brain structure abnormalities were excluded from the sample  
> T1 QA exclude: 23 participants with T1-weighted scans that failed visual QC were excluded from the sample  
> fMRI motion exclude: 179 participants with a mean relative RMS motion value > 0.2 during the resting state fMRI scan were excluded from the sample  
> Fluctuation amplitude outlier exclude: From the remaining 1052 participants, 19 individuals that had outlier (+- 4 SD from the mean) fluctuation amplitude data in more than 5% of glasser parcels were excluded from the sample   
<br>

### Model Fitting 

***GAM Functions***  
To characterize age-dependent changes in spontaneous activity fluctuations across the developing cortex as well as associations between fluctuation amplitude and environmental factors, generalized additive models were fit in each cortical region. GAM models and associated statistics, fitted values, smooths, and derivatives were quantified with the set of functions included in [/gam_models/GAM_functions.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/GAM_functions.R). This script includes the following functions:
-	*gam.fit.smooth*: A function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and save out statistics (smooth partial R squared and p-value) and derivative-based characteristics
-	*gam.smooth.predict*: A function to predict fitted values of a measure based on a GAM model and a prediction df
-	*gam.estimate.smooth*: A function to estimate zero-averaged gam smooth functions
-	*gam.posterior.smooths*: A function to simulate the posterior distribution from a fit GAM model, calculate smooths for individual posterior draws, and determine smooth max and min median values + 95% credible intervals
-	*gam.derivatives*: A function to compute derivatives for the smooth term from a main GAM model and for individual draws from the simulated posterior distribution; can return true model derivatives or posterior derivatives
-	*gam.fit.covariate*: A function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariate of interest + control covariates)) and save out statistics (partial R squared and p-value) for the first covariate

***Model Fitting: Age-Dependent Changes in Regional Fluctuation Amplitude***  
Developmental models were fit with [gam_models/fitGAMs_fluctuationamplitude_age.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/fitGAMs_fluctuationamplitude_age.R), which calls the functions described above. Age-focused GAMs were executed for Glasser 360 and Schaefer 400 atlases using the final study sample of N = 1033 generated during the Sample Construction process.  

***Model Fitting: Developmental Environment-Dependent Variation in Regional Fluctuation Amplitude***  
Models examining associations between fluctuation amplitude and environmental and family characteristics were run via [gam_models/fitGAMs_fluctuationamplitude_environment.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/fitGAMs_fluctuationamplitude_environment.R) on the N = 1033 study sample. 


### Data Interpretation and Visualization 

Model results were examined and studied within our hierarchical neurodevelopmental plasticity framework in the R markdown provided in [/developmental_effects/hierarchical_development.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd). This script does the following:

- Generates all manuscript Figures
- Examines the cortical distribution and significance of associations between regional fluctuation amplitude and age
- Quantifies associations between developmental changes in fluctuation amplitude and the T1w/T2w ratio
- Characterizes alignment of fluctuation amplitude age-related changes with the sensorimotor-association axis, including the magnitude and timing of change
- Calculates age smooths for 10 bins of the sensorimotor-association axis
- Performs a temporal sliding window analysis to uncover when developmental change in fluctuation amplitude is maximally aligned (and not aligned) with the sensorimotor-association axis
- Examines the significance and distribution of associations between regional fluctuation amplitude and environmental variability, and how environment effects are stratified by the S-A axis

### Sensitivity Analyses

The robustness of our developmental findings was evaluated in a series of sensitivity analyses. Sensitiviy analysis GAMs were fit with sensitivity_analyses/fitGAMs_sensitivityanalyses.R and results were examined in sensitivity_analyses/sensitivity_results.Rmd. The following sensitivity analyses were performed:

- *Low motion*: Findings were assessed in a low motion sample of N = 690 participants with a relative mean RMS < 0.075
- *Psychiatry exclusions*: Findings were assessed in a sample that excluded participants with current psychotropic medication use or a history of psychiatric hospitalization
- *CBF controlled*: All fluctuation amplitude models were re-fit while controlling for regional cerebral blood flow (quantified with ASL)
- *T2 controlled*: All fluctuation amplitude models were re-fit while controlling for regional BOLD signal level (T2*) during the fMRI scan
- *Normalization*: Developmental analyses were run with mean normalized fluctuation amplitude as the dependent variables
- *Atlas*: Results were reproduced using the Schaefer 400 atlas

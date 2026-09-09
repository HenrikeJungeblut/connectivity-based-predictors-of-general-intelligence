# Connectivity-based predictors of general intelligence

## File structure of the Code_package

```
Code_package/ 
│── Analysis_data/                                                                  # datasets created by the scripts in this code package 
│── atlas_labels/                                                                   # all values present in the atlas .nii image files for each subject and hemisphere
│   │── sub-010_restdata_labels_lh.txt
│   │── ...
│   └── sub-215_taskdata_labels_rh.txt
│── Data/                                                                                       
│   │── APM.json                                                                    # Advanced Progressive Matrices code book
│   │── APM.RData                                                                   # Advanced Progressive Matrices data frame
│   │── BIS_Complete_processed.csv                                                  # Berlin Intelligence Structure test data frame
│   │── BIS_Complete_processed.json                                                 # Berlin Intelligence Structure test code book
│   │── chanlocs.mat                                                                # for creating topoplot
│   │── df_full_EEG.Rda                                                             # behavioral data from EEG task
│   │── df_full_MRI.Rda                                                             # behavioral data from MRI task
│   └── subject_list.txt
│── erp/                                                                            # contains the ERP datasets
│   │── erp_ALS.csv                                                                 # manually picked P3 latencies for each participant and each task condition
│   │── switching_target_locked_sub-001.erp
│   │── ...
│   └── switching_target_locked_sub-215.erp
│── Figures/                                                                        # figures created by the scripts in this code package
│── MRI/
│   │── dwi/                                                                        # structural connectomes for each participant
│   │   │── sub-001_mean_FA_connectome.csv
│   │   │── ...
│   │   └── sub-215_mean_FA_connectome.csv
│   └── func/                                                                       # functional timeseries per participant and hemisphere of all Glasser parcels
│       │── rawfunc_timeseries_per_region/                                          # rawfunc: functional time series (task-evoked)
│       │    │── rawfunc_sub-010_lh_filtered_timeseries_per_glasser_region.txt          
│       │    │── ...
│       │    └── rawfunc_sub-215_rh_filtered_timeseries_per_glasser_region.txt
│       │── rest_timeseries_per_region/                                             # rest: resting-state time series
│       │    │── rest_sub-010_lh_filtered_timeseries_per_glasser_region.txt
│       │    │── ...
│       │    └── rest_sub-215_rh_filtered_timeseries_per_glasser_region.txt
│       └── timeseries_per_region/                                                  # res: residualized time series (background connectivity approach)
│            │── res_sub-010_lh_filtered_timeseries_per_glasser_region.txt
│            │── ...
│            └── res_sub-215_rh_filtered_timeseries_per_glasser_region.txt
│── Scripts/                                                                        # analysis scripts
│   │── EEG1_preprocessing.m                                                                   
│   │── ...                                                  
│   │── EEG4_plotting.m                                                                
│   │── MRI0_check_resampling_consequences.m  
│   │── ...                                                  
│   │── MRI8_SCFCcoupling_somatomotor.m                                                          
│   │── STATS0_Power_analysis.Rmd      
│   │── ...                                                 
│   │── STATS10_exploratory_analysis_CD.Rmd  
│   │── SUPP1_Descriptives_visual.Rmd      
│   │── ...                                                                                                  
│   └── SUPP6_SCFC_somatomotor_LC.Rmd
│── subject_connydegree_traces/                                                    # connectivity degree time series (based on PLI)                                                  
│   │── connydegree_traces_chan_Cz_pli_LM-R.csv                                    # for each channel and task condition
│   │── ...
│   └── connydegree_traces_chan_Y30_pli_OE-S.csv
│── Tools/                                                                         # necessary tools to run scripts
│   │── 2019_03_03_BCT/                                                            # Brain Connectivity Toolbox
│   │   └── ...
│   │── eeglab2025.0.0/                                                            # eeglab
│   │   └── ...
│   └── brewermap.m  
└──Code_package.Rproj                                                             # R-project to run all R-Markdown scripts in RStudio
```

## Instructions for running the scripts

- The scripts EEG1-EEG4 and MRI0-MRI8 are MATLAB scripts.
- The scripts STATS0-STATS10 and SUPP1-SUPP6 are R-Mardown files. Open the R-Project *Code_package.Rproj* in RStudio.

### Necessary software

The **Brain Connectivity Toolbox** and **eeglab** are provided in the *Tools/* directory.

The following R libraries need to be installed. We provide the version number of the packages used for this project's analysis.

| **Library** | **Version** |
| -------- | -----------|
| afex     | v.1.5-1    |
| dplyr    | v.1.2.1    |
| emmeans  | v.2.0.4    |
| ggplot2  | v.4.0.2    |
| ggpubr   | v.0.6.0    |
| here     | v.1.0.1    |
| lavaan   | v.0.6-21   |
| pheatmap | v.1.0.12   |
| psych    | v.2.4.6.26 |
| see      | v.0.12.0   |
| semPower | v.2.1.1    |
| tidyr    | v.1.3.1    |

## Order for running the scripts

### 1. EEG preprocessing and generation of EEG analysis data

1. EEG1_preprocessing.m **(NOT REPRODUCIBLE because raw EGG data is not provided in code package)**
2. EEG2_TF_CD_EKP.m **(NOT REPRODUCIBLE because raw EGG data is not provided in code package)**
3. EEG3_extract_CD_in_timewindows.m

### 2. Generation of MRI analysis data and MRI reliability

1. MRI0_check_resampling_consequences.m
2. MRI1_check_SC_matrices.m
3. MRI2_SCFCcoupling.m
4. MRI3_SCFCcoupling_Reliabilität.m
5. MRI4_SCFCcoupling_frontoparietal.m
6. MRI5_SCFCcoupling_frontoparietal_Reliabilität.m
7. MRI6_meanFC.m

### 3. Statistical analysis

1. STATS0_Power_analysis.Rmd
2. STATS1_Descriptives.Rmd
3. STATS2_Descriptives_rest.Rmd
4. STATS3_Descriptives_Beh.Rmd
5. STATS4_SCFC_whole_brain.Rmd
6. STATS5_SCFC_frontoparietal.Rmd
7. STATS6_FC_ANOVA.Rmd
8. STATS7_SCFC_frontoparietal_LC.Rmd
9. STATS8_CD.Rmd
10. STATS9_extract_latent_cors_CD.Rmd
11. EEG4_plotting.m
12. STATS10_exploratory_analysis_CD.Rmd

### 4. Supplementary analyses

1. MRI7_SCFCcoupling_visual.m
2. SUPP1_Descriptives_visual.Rmd
3. SUPP2_SCFC_visual.Rmd
4. SUPP3_SCFC_visual_LC.Rmd
5. MRI8_SCFCcoupling_somatomotor.m
6. SUPP4_Descriptives_somatomotor.Rmd
7. SUPP5_SCFC_somatomotor.Rmd
8. SUPP6_SCFC_somatomotor_LC.Rmd

## Code books for provided data sets

The following code books are provided as .json files in the *Data/* subdirectory:

1. APM.json
2. BIS_Complete_processed.json
3. df_full_EEG.json
4. df_full_MRI.json

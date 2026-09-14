# This script runs tedana (TE-dependent analysis) on multi-echo fMRI data 
# that was preprocessed by fMRIPrep. It takes the 4 individual echoes and 
# their specific echo times (TEs) to automatically distinguish between BOLD 
# signal and non-BOLD noise, creating an optimally combined, denoised image.
for sub in 001
do

cd /path/to/derivative/sub-${sub}/fmriprep_switch/ses-01/func/ #adjust path

tedana -d sub-${sub}_ses-01_task-switch_acq-mb3me4_echo-1_desc-preproc_bold.nii.gz \
 sub-${sub}_ses-01_task-switch_acq-mb3me4_echo-2_desc-preproc_bold.nii.gz \
 sub-${sub}_ses-01_task-switch_acq-mb3me4_echo-3_desc-preproc_bold.nii.gz \
 sub-${sub}_ses-01_task-switch_acq-mb3me4_echo-4_desc-preproc_bold.nii.gz \
 -e 12 28.24 44.48 60.72 \
#adjust path
 --out-dir /path/to/derivative/sub-${sub}/tedana_switch/ \
 --mask sub-${sub}_ses-01_task-switch_acq-mb3me4_desc-brain_mask.nii.gz \
 --overwrite

done

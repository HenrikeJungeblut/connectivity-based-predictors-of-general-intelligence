# This script prepares DWI data for susceptibility distortion correction using 
# FSL's topup. It extracts the b0 volumes from AP and PA phase-encoded scans, 
# merges them, calculates the distortion field with topup, and averages the 
# corrected b0 images to create a single high-fidelity reference image.
cat /path/to/subject_names.txt | tr -d '\r' | while read sub; do

# create output dir if it doesnt exist
# mkdir /path/to/derivative/sub-${sub}
# mkdir /path/to/derivative/sub-${sub}/dwi

# fslroi cuts images into pieces --> we cut our images into pieces in terms of the volumes --> takes the first volume (b0 image)
fslroi /path/to/raw_data/${sub}/dwi/${sub}_dwi_AP.nii.gz /path/to/derivative/${sub}/dwi/nodif.nii.gz 0 1

fslroi /path/to/raw_data/${sub}/dwi/${sub}_dwi_PA.nii.gz /path/to/derivative/${sub}/dwi/nodif_PA.nii.gz 0 1

fslmerge -t /path/to/derivative/${sub}/dwi/AP_PA_b0.nii.gz /path/to/derivative/${sub}/dwi/nodif.nii.gz \
/path/to/derivatives/${sub}/dwi/nodif_PA.nii.gz

topup --imain=/path/to/derivative/${sub}/dwi/AP_PA_b0 --datain=/path/to/derivative/2_scripts_DWI/acq_param.txt \
--config=b02b0.cnf --out=/path/to/derivative/${sub}/dwi/topup_AP_PA_b0 \
--fout=/path/to/derivative/${sub}/dwi/topup_AP_PA_b0_fout \
--iout=/path/to/derivative/${sub}/dwi/topup_AP_PA_b0_iout

# topup_warp is the distortion-corrected image

# last step is used to create a brain mask by computing the average image of the corrected b0
fslmaths /path/to/derivative/${sub}/dwi/topup_AP_PA_b0_iout.nii.gz \
-Tmean /path/to/derivative/${sub}/dwi/hifi_nodif

echo 'Done with Subject '$sub
done

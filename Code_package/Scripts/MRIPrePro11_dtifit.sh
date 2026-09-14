# This script prepares files for tensor fitting and runs FSL's dtifit. It 
# extracts and skull-strips the eddy-corrected b0 image, then calculates 
# DTI scalar maps (like FA and MD) using the eddy-corrected diffusion data 
# and the updated, eddy-rotated gradient vectors (bvecs).
cat /path/to/subject_names.txt | tr -d '\r' | while read sub; do

echo "working on $sub"
# mask needs to be renamed to nodif_brain_mask so that dtifit will work
cp /path/to/derivative/${sub}/dwi/topup_AP_PA_b0_iout_mean_brain_mask.nii.gz \
/path/to/derivative/${sub}/dwi/nodif_brain_mask.nii.gz

fslroi /path/to/derivative/${sub}/dwi/eddy_corrected_data.nii.gz \
/path/to/derivative/${sub}/dwi/eddy_corrected_data_b0.nii.gz 0 1

fslmaths /path/to/derivative/${sub}/dwi/eddy_corrected_data_b0.nii.gz \
-mul /path/to/derivative/${sub}/dwi/nodif_brain_mask.nii.gz \
/path/to/derivative/${sub}/dwi/nodif_brain.nii.gz

dtifit --data=/path/to/derivative/${sub}/dwi/eddy_corrected_data.nii.gz \
--out=/path/to/derivative/${sub}/dwi/$sub \
--mask=/path/to/derivative/${sub}/dwi/nodif_brain_mask.nii.gz \
--bvecs=/path/to/derivative/${sub}/dwi/eddy_corrected_data.eddy_rotated_bvecs \
--bvals=/path/to/derivative/${sub}/dwi/bvals

echo 'Done SUB-'$sub
done


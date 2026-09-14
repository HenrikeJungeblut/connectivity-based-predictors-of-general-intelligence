# This script performs brain extraction (skull-stripping) on the distortion-
# corrected b0 image using FSL's bet to generate a binary brain mask. It also 
# copies and renames the gradient tables (bvals/bvecs) to the standard naming 
# convention required by downstream FSL tools like eddy.

cat /path/to/derivative/subject_names.txt | tr -d '\r' | while read sub; do

# -o generates brain surface outline overlaid onto original image
# -m creates binary brain mask
# choose coordinates from voxel location
# 0.3 is somewhat smaller than the default (0.5), but the default often cuts away frontal parts of the brain
# value needs to be found for every subject
bet /path/to/derivative/${sub}/dwi/hifi_nodif.nii.gz /path/to/derivative/${sub}/dwi/hifi_nodif_brain.nii.gz -f 0.3 -c 70 70 42 -o -m

# renaming bvecs and bvals files so that FSL can work with them 
cp /path/to/raw_data/${sub}/dwi/${sub}_dwi_AP.bvec /path/to/derivative/${sub}/dwi/bvecs

cp /path/to/raw_data/${sub}/dwi/${sub}_dwi_AP.bval /path/to/derivative/${sub}/dwi/bvals

done



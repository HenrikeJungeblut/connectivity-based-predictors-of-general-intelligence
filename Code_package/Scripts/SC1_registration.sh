# This script converts the surface-based Glasser (HCP-MMP1) atlas into 3D volumes 
# in structural space, and then warps those volumes into native diffusion (DTI) 
# space. It uses FreeSurfer to map the surface labels to the T1 volume, and 
# MRtrix3 to apply the previously calculated T1-to-DTI transformation (MRIPrePro3).

cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
  echo "registering $SUBJECT"

mri_label2vol \
--subject $SUBJECT \
--hemi lh \
--annot HCP-MMP1 \
--temp /path/to/derivative/freesurfer/$SUBJECT/mri/brain.mgz \
--o /path/to/derivative/freesurfer/$SUBJECT/mri/lh.HCP-MMP1.nii.gz \
--identity

mri_label2vol \
--subject $SUBJECT \
--hemi rh \
--annot HCP-MMP1 \
--temp /path/to/derivative/freesurfer/$SUBJECT/mri/brain.mgz \
--o /path/to/derivative/freesurfer/$SUBJECT/mri/rh.HCP-MMP1.nii.gz \
--identity

# Define base path
BASE_PATH="/path/to"
DERIV="$BASE_PATH/derivative/freesurfer"
RAW="$BASE_PATH/raw_data/$SUBJECT/anat"

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrtransform \
  "$DERIV/${SUBJECT}/mri/lh.HCP-MMP1.nii.gz" \
  -linear "$RAW/${SUBJECT}_T1w2dti_mrtrixformat.txt" \
  "$DERIV/${SUBJECT}/mri/lh.HCP-MMP1_dti.nii.gz"    \
  -inverse \
  -force

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrtransform \
  "$DERIV/${SUBJECT}/mri/rh.HCP-MMP1.nii.gz" \
  -linear "$RAW/${SUBJECT}_T1w2dti_mrtrixformat.txt" \
  "$DERIV/${SUBJECT}/mri/rh.HCP-MMP1_dti.nii.gz"    \
  -inverse \
  -force
done

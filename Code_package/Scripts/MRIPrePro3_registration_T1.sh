# This script registers the structural T1w brain image to the native diffusion 
# (DTI) space. It uses FSL's flirt to compute the linear transformation, converts 
# the matrix to MRtrix format via Docker, and then applies the inverse of that 
# transformation to map the T1w brain into the DTI space using mrtransform.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
  echo "registering $SUBJECT"

# register T1w image to dti 
/usr/local/fsl/fslv6/bin/flirt -in /path/to/derivative/$SUBJECT/dwi/eddy_corrected_data_biascorrected.nii.gz \
-ref /path/to/raw_data/$SUBJECT/anat/${SUBJECT}_T1w_brain.nii.gz \
-omat /path/to/raw_data/$SUBJECT/anat/${SUBJECT}_T1w2dti.mat \
-bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear

# Define base path
BASE_PATH="/path/to/"
DERIV="$BASE_PATH/derivative/$SUBJECT/dwi"
RAW="$BASE_PATH/raw_data/$SUBJECT/anat"


docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  transformconvert \
  "$RAW/${SUBJECT}_T1w2dti.mat" \
  "$DERIV/eddy_corrected_data_biascorrected.nii.gz" \
  "$RAW/${SUBJECT}_T1w.nii.gz" \
  flirt_import "$RAW/${SUBJECT}_T1w2dti_mrtrixformat.txt" \
  -force

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrtransform \
  "$RAW/${SUBJECT}_T1w_brain.nii.gz" \
  -linear "$RAW/${SUBJECT}_T1w2dti_mrtrixformat.txt" \
  "$RAW/${SUBJECT}_T1w2dti.nii.gz"    \
  -inverse \
  -force

done



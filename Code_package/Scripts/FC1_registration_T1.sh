# This script registers the structural T1w brain image to the native functional 
# (fMRI) space. It uses FSL's flirt to compute the linear transformation from the 
# tedana-denoised functional data to the T1w, converts the matrix to MRtrix format, 
# and applies the inverse transformation to map the T1w brain into functional space.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
  echo "registering $SUBJECT"

# register T1w image to func 
/usr/local/fsl/fslv6/bin/flirt -in /path/to/derivative/$SUBJECT/tedana_switch/desc-denoised_bold.nii.gz \
-ref /path/to/raw_data/$SUBJECT/anat/${SUBJECT}_T1w_brain.nii.gz \
-omat /path/to/raw_data/$SUBJECT/anat/${SUBJECT}_T1w2func.mat \
-bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear

# Define base path
BASE_PATH="/path/to"
DERIV="$BASE_PATH/derivative/$SUBJECT/tedana_switch"
RAW="$BASE_PATH/raw_data/$SUBJECT/anat"


docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  transformconvert \
  "$RAW/${SUBJECT}_T1w2func.mat" \
  "$DERIV/desc-denoised_bold.nii.gz" \
  "$RAW/${SUBJECT}_T1w.nii.gz" \
  flirt_import "$RAW/${SUBJECT}_T1w2func_mrtrixformat.txt" \
  -force

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrtransform \
  "$RAW/${SUBJECT}_T1w_brain.nii.gz" \
  -linear "$RAW/${SUBJECT}_T1w2func_mrtrixformat.txt" \
  "$RAW/${SUBJECT}_T1w2func.nii.gz"    \
  -inverse \
  -force

done



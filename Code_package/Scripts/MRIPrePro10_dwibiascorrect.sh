# This script loops through subjects to perform B1 bias field correction on 
# eddy-corrected diffusion data. It uses the ANTs N4 algorithm via MRtrix3's 
# dwibiascorrect, run safely inside a Docker container with local volume mounts 
# and a dedicated temporary scratch directory that is cleaned up afterward.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
DERIV="$BASE_PATH/derivative/$SUBJECT/dwi"
RAW="$BASE_PATH/raw_data/$SUBJECT/dwi"
SCRATCH_DIR="$DERIV/tmp"

# Make sure scratch dir exists
mkdir -p "$SCRATCH_DIR"

docker run --rm -i \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 dwibiascorrect ants \
  "$DERIV/eddy_corrected_data.nii.gz" \
  "$DERIV/eddy_corrected_data_biascorrected.nii.gz" \
  -fslgrad "$RAW/${SUBJECT}_dwi_AP.bvec" "$RAW/${SUBJECT}_dwi_AP.bval" \
  -mask "$DERIV/hifi_nodif_brain_mask.nii.gz" \
  -scratch "$SCRATCH_DIR" \
  -force

rm -rf "$SCRATCH_DIR"

done


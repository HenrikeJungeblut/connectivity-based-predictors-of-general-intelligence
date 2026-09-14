# This script converts the fully preprocessed DWI data into MRtrix format (.mif) 
# and estimates the multi-tissue response functions (white matter, grey matter, 
# and CSF) using the Dhollander algorithm, which is required for downstream 
# Constrained Spherical Deconvolution (CSD).
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
DERIV="$BASE_PATH/derivative/$SUBJECT/dwi"
RAW="$BASE_PATH/raw_data/$SUBJECT/dwi"
SCRATCH_DIR="$DERIV/tmp"

mkdir -p "$SCRATCH_DIR"

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrconvert \
  -fslgrad "$DERIV/eddy_corrected_data.eddy_rotated_bvecs" "$DERIV/bvals" \
  "$DERIV/eddy_corrected_data_biascorrected.nii.gz" \
  "$DERIV/dwi.mif" \
  -force


docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -w "$DERIV" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 dwi2response dhollander \
  "$DERIV/dwi.mif" \
  "$DERIV/out_sfwm.txt" "$DERIV/out_gm.txt" "$DERIV/out_csf.txt" \
  -mask "$DERIV/nodif_brain_mask_dilated.mif" \
  -fslgrad "$DERIV/eddy_corrected_data.eddy_rotated_bvecs" "$DERIV/bvals" \
  -force

rm -rf "$SCRATCH_DIR"

done

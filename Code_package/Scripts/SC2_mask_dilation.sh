# This script uses MRtrix3's 'maskfilter' via Docker to dilate (expand) the 
# boundaries of the diffusion brain mask. 
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
DERIV="$BASE_PATH/derivative/$SUBJECT/dwi"
RAW="$BASE_PATH/raw_data/$SUBJECT/anat"

docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 maskfilter \
  "$DERIV/nodif_brain_mask.nii.gz" \
  dilate \
  "$DERIV/nodif_brain_mask_dilated.mif" 


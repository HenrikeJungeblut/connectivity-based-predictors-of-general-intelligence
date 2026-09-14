# This script performs Multi-Shell Multi-Tissue Constrained Spherical 
# Deconvolution (MSMT-CSD) using MRtrix3. It takes the multi-tissue response 
# functions and the preprocessed DWI data to estimate the Fiber Orientation 
# Distributions (FODs) for white matter, grey matter, and CSF.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
DERIV="$BASE_PATH/derivative/$SUBJECT/dwi"
RAW="$BASE_PATH/raw_data/$SUBJECT/dwi"
SCRATCH_DIR="$DERIV/tmp"

#mkdir -p "$SCRATCH_DIR"

docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -w "$DERIV" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 dwi2fod msmt_csd \
  "$DERIV/dwi.mif" \
  "$DERIV/out_sfwm.txt" "$DERIV/wmfod.mif" "$DERIV/out_gm.txt" "$DERIV/gmfod.mif" "$DERIV/out_csf.txt" "$DERIV/csffod.mif" 

done

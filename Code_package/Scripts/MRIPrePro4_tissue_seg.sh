# This script generates a 5-Tissue-Type (5TT) segmentation image from the 
# DTI-aligned T1w image, which is required for Anatomically-Constrained 
# Tractography (ACT). It uses MRtrix3's '5ttgen fsl' via Docker to create the 
# segmentation, and then uses 'mrconvert' to save a NIfTI copy of the result.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
RAW="$BASE_PATH/raw_data/$SUBJECT/anat"
DERIV="$BASE_PATH/derivative/$SUBJECT/anat"
SCRATCH_DIR="$DERIV/tmp"

mkdir -p "$DERIV"
mkdir -p "$SCRATCH_DIR"

docker run --rm -it \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  5ttgen fsl \
  "$RAW/${SUBJECT}_T1w2dti.nii.gz" \
  "$DERIV/5TT.mif" \
  -scratch "$SCRATCH_DIR" \
  -premask \
  -noclean \
  -force

rm -rf "$SCRATCH_DIR"

# convert to visualize in viewer for quality control
docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrconvert \
  "$DERIV/5TT.mif" \
  "$DERIV/5TT.nii.gz" 

done

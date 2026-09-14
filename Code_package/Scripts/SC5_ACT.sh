# This script performs whole-brain probabilistic tractography using MRtrix3. 
# It generates 30 million streamlines based on the white matter FODs, utilizing 
# Anatomically-Constrained Tractography (ACT), dynamic seeding, and crops the 
# endpoints exactly at the grey matter-white matter interface.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
DWI="$BASE_PATH/derivative/$SUBJECT/dwi"
ANAT="$BASE_PATH/derivative/$SUBJECT/anat"
FS="$BASE_PATH/derivative/freesurfer/$SUBJECT/mri"
SCRATCH_DIR="$DERIV/tmp"

#mkdir -p "$SCRATCH_DIR"

docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -w "$DERIV" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 tckgen \
  "$DWI/wmfod.mif" \
  "$DWI/tracks.tck" \
  -act "$ANAT/5TT.mif" \
  -backtrack \
  -seed_dynamic "$DWI/wmfod.mif" \
  -crop_at_gmwmi \
  -select 30000000 \
  -force

done

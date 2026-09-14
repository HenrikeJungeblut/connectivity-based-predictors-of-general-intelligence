# This script filters the whole-brain tractogram using MRtrix3's SIFT algorithm. 
# It selectively removes streamlines so that the density of the final tractogram 
# accurately matches the underlying apparent fiber density estimated from the 
# white matter FODs, which is a crucial step for quantitative connectomics.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
DWI="$BASE_PATH/derivative/$SUBJECT/dwi"
ANAT="$BASE_PATH/derivative/$SUBJECT/anat"

#mkdir -p "$SCRATCH_DIR"

docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -w "$DERIV" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 tcksift \
  -act "$ANAT/5TT.mif" \
  "$DWI/tracks.tck" \
  "$DWI/wmfod.mif" \
  "$DWI/filtered_tracks.tck" \
  -force 

done

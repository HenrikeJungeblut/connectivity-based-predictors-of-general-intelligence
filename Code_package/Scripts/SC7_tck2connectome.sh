# This script constructs a structural connectivity matrix weighted by Fractional 
# Anisotropy (FA). It first combines the left and right hemisphere Glasser atlases 
# (adding a 180-label offset to the right hemisphere for unique IDs). It then 
# samples the mean FA along every streamline and uses tck2connectome to build 
# a symmetric, FA-weighted 360x360 connectome matrix.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do

# Define base path
BASE_PATH="/path/to"
DWI="$BASE_PATH/derivative/$SUBJECT/dwi"
ANAT="$BASE_PATH/derivative/$SUBJECT/anat"
OUT="$BASE_PATH/derivative/scripts/4_SC-FC-coupling/DWI/mean_FA_connectome"
FS="$BASE_PATH/derivative/freesurfer/$SUBJECT/mri"
SCRATCH_DIR="$DERIV/tmp"

#mkdir -p "$SCRATCH_DIR"

## prepare parcellation
# convert to mif
docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrconvert \
  "$FS/lh.HCP-MMP1_dti.nii.gz" \
  "$FS/lh.HCP-MMP1_dti.mif" \
  -force

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrconvert \
  "$FS/rh.HCP-MMP1_dti.nii.gz" \
  "$FS/rh.HCP-MMP1_dti.mif" \
  -force 

#add offset and combine both hemispheres
docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrcalc \
  "$FS/rh.HCP-MMP1_dti.mif" \
  0 -gt \
  180 -mul \
  "$FS/rh.HCP-MMP1_dti.mif" -add \
  "$FS/rh.HCP-MMP1_dti_offset.mif" \
  -force

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrcalc \
  "$FS/lh.HCP-MMP1_dti.mif" "$FS/rh.HCP-MMP1_dti_offset.mif" -add "$FS/HCP-MMP1_dti.mif" \
  -force

docker run --rm -it \
  -u $(id -u ${USER}):2018 \
  -v "$BASE_PATH":"$BASE_PATH" \
  mrtrix3/mrtrix3:latest \
  mrconvert \
  "$DWI/${SUBJECT}_FA.nii.gz" \
  "$DWI/${SUBJECT}_FA.mif" \
  -force 

docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -w "$DERIV" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 tcksample \
  "$DWI/filtered_tracks.tck" \
  "$DWI/${SUBJECT}_FA.mif" \
  "$DWI/mean_FA_per_streamline.csv" \
  -stat_tck mean \
  -force

docker run --rm -it \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  -v "$BASE_PATH":"$BASE_PATH" \
  -w "$DERIV" \
  -e MRTRIX_TMPFILE_DIR="$SCRATCH_DIR" \
  mrtrix3/mrtrix3 tck2connectome \
  "$DWI/filtered_tracks.tck" \
  "$FS/HCP-MMP1_dti.mif" \
  "$OUT/${SUBJECT}_mean_FA_connectome.csv" \
  -scale_file "$DWI/mean_FA_per_streamline.csv" \
  -stat_edge mean \
  -symmetric -zero_diagonal \
  -force

done

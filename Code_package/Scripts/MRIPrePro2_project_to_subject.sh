# This script loops through a list of subject IDs and uses FreeSurfer's 
# mri_surf2surf to project the Glasser (HCP-MMP1) cortical atlas from the 
# standard fsaverage template onto each individual's native brain surface.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
  echo "Projecting Glasser atlas to $SUBJECT"

  # Left hemisphere
  mri_surf2surf --srcsubject fsaverage --trgsubject $SUBJECT \
    --hemi lh --sval-annot /path/to/glasser_atlas/lh.HCP-MMP1.annot \
    --tval $SUBJECTS_DIR/$SUBJECT/label/lh.HCP-MMP1.annot

  # Right hemisphere
  mri_surf2surf --srcsubject fsaverage --trgsubject $SUBJECT \
    --hemi rh --sval-annot /path/to/glasser_atlas/rh.HCP-MMP1.annot \
    --tval $SUBJECTS_DIR/$SUBJECT/label/rh.HCP-MMP1.annot
done



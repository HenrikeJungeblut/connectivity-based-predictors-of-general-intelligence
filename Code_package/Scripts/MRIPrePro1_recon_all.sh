# cortical reconstruction and volumetric segmentation of the T1-weighted data performed using Freesurfer's recon-all
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
  recon-all -s $SUBJECT -i /path/to/raw_data/$SUBJECT/anat/${SUBJECT}_T1w.nii.gz -all
done


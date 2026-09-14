# This script extracts a list of all unique region IDs (labels) that are actually 
# present in the subject's functional space atlases. It converts the NIfTI volumes 
# to text using fsl2ascii, then sorts the voxel values to identify exactly which 
# Glasser regions survived the parcellation and resampling steps.
cat /path/to/subject_names.txt | tr -d '\r' | while read SUBJECT; do
# task
## left 
fsl2ascii /path/to/derivative/freesurfer/$SUBJECT/mri/HCP-MMP1_task_cortexonly_lh_RS.nii.gz /path/to/atlas_labels/${SUBJECT}_taskdata_labels_lh.txt
cat /path/to/atlas_labels/${SUBJECT}_taskdata_labels_lh* | tr ' ' '\n' | grep -v '^$' | sort -u  > /path/to/atlas_labels/${SUBJECT}_taskdata_labels_lh.txt
## right
fsl2ascii /path/to/derivative/freesurfer/$SUBJECT/mri/HCP-MMP1_task_cortexonly_rh.nii.gz /path/to/atlas_labels/${SUBJECT}_taskdata_labels_rh.txt
cat /path/to/atlas_labels/${SUBJECT}_taskdata_labels_rh* | tr ' ' '\n' | grep -v '^$' | sort -u  > /path/to/atlas_labels/${SUBJECT}_taskdata_labels_rh.txt

# rest
## left
fsl2ascii /path/to/derivative/freesurfer/$SUBJECT/mri/HCP-MMP1_rest_cortexonly_lh.nii.gz /path/to/atlas_labels/${SUBJECT}_restdata_labels_lh.txt
cat /path/to/atlas_labels/${SUBJECT}_restdata_labels_lh* | tr ' ' '\n' | grep -v '^$' | sort -u  > /path/to/atlas_labels/${SUBJECT}_restdata_labels_lh.txt
## right
fsl2ascii /path/to/derivative/freesurfer/$SUBJECT/mri/HCP-MMP1_rest_cortexonly_rh.nii.gz /path/to/atlas_labels/${SUBJECT}_restdata_labels_rh.txt
cat /path/to/atlas_labels/${SUBJECT}_restdata_labels_rh* | tr ' ' '\n' | grep -v '^$' | sort -u  > /path/to/atlas_labels/${SUBJECT}_restdata_labels_rh.txt
done

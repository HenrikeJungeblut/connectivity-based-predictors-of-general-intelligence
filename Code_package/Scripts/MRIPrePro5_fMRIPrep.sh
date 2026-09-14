# This script runs the fMRIPrep preprocessing pipeline via Docker for a specific 
# task ("switch"). It is configured to output individual echoes for multi-echo 
# data and specifically bypasses BIDS validation and FreeSurfer surface 
# reconstruction (--fs-no-reconall) to speed up the processing.

list="sub-001"

#########################################################
#Enter path to original bids base directory
BidsBase="/path/to"


for d in $list ; do
# Enter docker command
	docker run -ti --rm -u $(id -u ${USER}):$(id -g ${USER}) \
	-v $BidsBase:/base \
	nipreps/fmriprep:latest \
	/base/raw_data \
	/base/derivative \
	participant \
	--participant-label $d \
	--fs-license-file /base/license.txt \
	--me-output-echos \
	--no-submm-recon \
	--fs-no-reconall \
	--skip-bids-validation \
	-t switch

done


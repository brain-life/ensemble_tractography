#!/bin/bash
## Cadidate fascicles generation script.
## 
## This shell script uses mrtrix/0.2.12 to run a single tractography method using probabilistic
## tractography based on the constrained spherical deconvolution. 
##
## Brent McPherson and Franco Pestilli Indiana University 2016

## load necessary modules on the Karst cluster environment
#module switch mrtrix/0.2.12

## The script requires a single input tha tis the folder name for the subject that needs to be processed
## The corrent version of the script handles only data on the local Indiana University Cluster Systems
## Under the project lifebid. There are two current data sets there the HCP and STN96
SUBJ=$1

## set lmax
i_lmax=8

## set tracking
i_tracktype=SD_PROB

## Set paths to diffusion data directories
## DWIFILENAME=run01_fliprot_aligned_trilin
## TOPDIR=/N/dc2/projects/lifebid/2t1/predator/$SUBJ

## HCP dataset
DWIFILENAME=6_DWI
TOPDIR=/N/dc2/projects/lifebid/glue/subjects/$SUBJ

## STN96 dataset
ANATDIR=$TOPDIR/anat
OUTDIR=$TOPDIR/fibers

mkdir -v $OUTDIR

## Number of fibers requested and max number attempted to hit the number.
NUMFIBERS=500000
MAXNUMFIBERSATTEMPTED=1000000

echo 
echo Performing preprocessing of data before starting tracking...
echo 

## convert wm mask
mrconvert $ANATDIR/wm_mask.nii.gz $OUTDIR/${DWIFILENAME}_wm.mif

## convert dwi's
bash /N/dc2/projects/lifebid/glue/bin/convert_nifti_mrtrix.sh $SUBJ

## convert dwi's 
#mrconvert $TOPDIR/diffusion/${DWIFILENAME}.nii.gz $OUTDIR/${DWIFILENAME}.mif

## make mask from DWI data
#average $OUTDIR/${DWIFILENAME}.mif -axis 3 - | threshold - - | median3D - - | median3D - $OUTDIR/${DWIFILENAME}_brainmask.mif

## make bet mask because this sucks
#bet $TOPDIR/diffusion/${DWIFILENAME}.nii.gz $TOPDIR/diffusion/${DWIFILENAME}_bet -R -m -f 0.3
mrconvert $TOPDIR/diffusion/${DWIFILENAME}_bet_mask.nii.gz $OUTDIR/${DWIFILENAME}_brainmask.mif

## make gradient
#cp $TOPDIR/diffusion/$DWIFILENAME.b $OUTDIR
#paste $TOPDIR/diffusion/$DWIFILENAME.bvec $TOPDIR/diffusion/$DWIFILENAME.bval > $OUTDIR/tmp.b
#sed "s/\,/\\t/g" $OUTDIR/tmp.b > $OUTDIR/$DWIFILENAME.b
#rm $OUTDIR/tmp.b

#cat $OUTDIR/$DWIFILENAME.b | tr ',' ' ' > $OUTDIR/$DWIFILENAME.b

## fit tensors
dwi2tensor $OUTDIR/${DWIFILENAME}.mif -grad $OUTDIR/$DWIFILENAME.b $OUTDIR/${DWIFILENAME}_dt.mif 

## create FA image
tensor2FA $OUTDIR/${DWIFILENAME}_dt.mif - | mrmult - $OUTDIR/${DWIFILENAME}_brainmask.mif $OUTDIR/${DWIFILENAME}_fa.mif

##
## CREATE TEMPORARY FAKE WM MASK
##

## create a whole brain mask from diffusion data
#threshold $OUTDIR/${DWIFILENAME}_fa.mif - -abs 0.2 | erode - $OUTDIR/${DWIFILENAME}_wm.mif

## create eigenvector map
tensor2vector $OUTDIR/${DWIFILENAME}_dt.mif - | mrmult - $OUTDIR/${DWIFILENAME}_fa.mif $OUTDIR/${DWIFILENAME}_ev.mif

## # Estimate deconvolution kernel: Estimate the kernel for deconvolution, using voxels with highest FA
## erodes brainmask - removes extreme artifacts (w/ high FA), creates FA image, AND single fiber mask 
erode $OUTDIR/${DWIFILENAME}_brainmask.mif -npass 3 - | mrmult $OUTDIR/${DWIFILENAME}_fa.mif - - | threshold - -abs 0.7 $OUTDIR/${DWIFILENAME}_sf.mif

## estimates response function
estimate_response $OUTDIR/${DWIFILENAME}.mif $OUTDIR/${DWIFILENAME}_sf.mif -lmax $i_lmax -grad $OUTDIR/$DWIFILENAME.b $OUTDIR/${DWIFILENAME}_response.txt
## # End estimation of deconvolution kernel

## Perform CSD in each white matter voxel
csdeconv $OUTDIR/${DWIFILENAME}.mif -grad $OUTDIR/$DWIFILENAME.b $OUTDIR/${DWIFILENAME}_response.txt -lmax $i_lmax -mask $OUTDIR/${DWIFILENAME}_brainmask.mif $OUTDIR/${DWIFILENAME}_lmax${i_lmax}.mif

# for i_lmax in 2 4 6 8; do
#     csdeconv $OUTDIR/${DWIFILENAME}.mif -grad $OUTDIR/$DWIFILENAME.b $OUTDIR/${DWIFILENAME}_response.txt -lmax $i_lmax -mask $OUTDIR/${DWIFILENAME}_brainmask.mif $OUTDIR/${DWIFILENAME}_lmax${i_lmax}.mif
# ## echo DONE Lmax=$i_lmax 
# done 

echo DONE performing preprocessing of data before starting tracking...
echo START tracking...

## track 1 connectome
streamtrack $i_tracktype $OUTDIR/${DWIFILENAME}_lmax${i_lmax}.mif $OUTDIR/${DWIFILENAME}_csd_lmax${i_lmax}_wm_${i_tracktype}-NUM${i_track}-$NUMFIBERS.tck -seed $OUTDIR/${DWIFILENAME}_wm.mif -mask $OUTDIR/${DWIFILENAME}_wm.mif -grad $OUTDIR/$DWIFILENAME.b -number $NUMFIBERS -maxnum $MAXNUMFIBERSATTEMPTED

# ##echo tracking Deterministic Tensorbased
# for i_track in 01 02 03 04 05 06 07 08 09 10; do
# streamtrack DT_STREAM $OUTDIR/${DWIFILENAME}.mif $OUTDIR/${DWIFILENAME}_wm_tensor-NUM${i_track}-$NUMFIBERS.tck -seed $OUTDIR/${DWIFILENAME}_wm.mif -mask $OUTDIR/${DWIFILENAME}_wm.mif -grad $OUTDIR/${DWIFILENAME}.b -number $NUMFIBERS -maxnum $MAXNUMFIBERSATTEMPTED
# done

# for i_track in 01 02 03 04 05 06 07 08 09 10; do
# ## loop over tracking and lmax
# for i_tracktype in SD_STREAM SD_PROB; do
# ##
# ##echo Tracking $i_tracktype Deterministic=1 Probabilistic=2 CSD-based
# ##
#     for i_lmax in 2 4 6 8; do
# 	##echo Tracking CSD-based Lmax=$i_lmax
# 	streamtrack $i_tracktype $OUTDIR/${DWIFILENAME}_lmax${i_lmax}.mif   $OUTDIR/${DWIFILENAME}_csd_lmax${i_lmax}_wm_${i_tracktype}-NUM${i_track}-$NUMFIBERS.tck -seed $OUTDIR/${DWIFILENAME}_wm.mif  -mask $OUTDIR/${DWIFILENAME}_wm.mif  -grad $OUTDIR/$DWIFILENAME.b -number $NUMFIBERS -maxnum $MAXNUMFIBERSATTEMPTED
#     done
# ##
# ##echo DONE Tracking $i_tracktype Deterministic=1 Probabilistic=2 CSD-based
# ##
# done
# done

echo DONE tracking. 

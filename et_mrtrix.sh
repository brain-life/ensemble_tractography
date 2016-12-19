#!/bin/bash
## Ensemble tractography.
## 
## This shell script uses mrtrix/0.2.12 to run a series of tractography methods using both probabilistic
## and deterministic tractography based on the tensor model or on constrained spherical deconvolution. 
##
## Brent McPherson and Franco Pestilli Indiana University 2016

## Make sure mrtrix is installed.
## 
## Locally we load mrtrix as a module on our clusters. 
module load mrtrix/0.2.12

## The script requires a single input that is the folder name for the subject that needs to be processed.
##
SUBJ=$1

## Set paths to diffusion data directories, the anatomy, and the outputdirectory where all files will be saved.
DWIFILENAME=<this_is_the_dwi_file_name>
TOPDIR=/full/path/to/subject/dwi/data/$SUBJ
ANATDIR=$TOPDIR/anatomy/
OUTDIR=$TOPDIR/ensemble_tractograms

mkdir -v $OUTDIR

## Number of fascicles requested and max number attempted to hit the number.
NUMFIBERS=500000
MAXNUMFIBERSATTEMPTED=1000000

##
echo 
echo Performing preprocessing of data before starting tracking...
echo 
##

## A good white matter mask is necessary for tracking. This is a NIFTI file with 1's in all voxels classified as WM.
## The mask is generally stored in the anatomy folder. Hereafter we assume the mask was properly created and 
## convert it into a mif file
mrconvert $ANATDIR/wm_mask.nii.gz $OUTDIR/${DWIFILENAME}_wm.mif

## We also convert dwi's nifiti files into mif files. 
mrconvert $TOPDIR/diffusion_data/$DWIFILENAME.nii.gz $OUTDIR/${DWIFILENAME}_dwi.mif

## We next want to estimate a response function to use for constrained-spherical deconvolution.
##
## To do so we find locations of high-FA (fractional anisotropy) in the brain and estimate the 
## CSD response in those voxels assuming that they contain a "single fiber" (e.g., the corpus callusom).
##
## We then make rough brain-mask from the dwi data.
average $OUTDIR/${DWIFILENAME}_dwi.mif -axis 3 - | threshold - - | median3D - - | median3D - $OUTDIR/${DWIFILENAME}_brainmask.mif

## Fit the tensor model
dwi2tensor $OUTDIR/${DWIFILENAME}_dwi.mif -grad $OUTDIR/$DWIFILENAME.b $OUTDIR/${DWIFILENAME}_dt.mif 

## Compute FA only in voxels within the brain.
tensor2FA $OUTDIR/${DWIFILENAME}_dt.mif - | mrmult - $OUTDIR/${DWIFILENAME}_brainmask.mif $OUTDIR/${DWIFILENAME}_fa.mif

## Create and principak diffusion direction (RGB) map
tensor2vector $OUTDIR/${DWIFILENAME}_dt.mif - | mrmult - $OUTDIR/${DWIFILENAME}_fa.mif $OUTDIR/${DWIFILENAME}_ev.mif

## erodes brainmask to find voxels with high FA. 
erode $OUTDIR/${DWIFILENAME}_brainmask.mif -npass 3 - | mrmult $OUTDIR/${DWIFILENAME}_fa.mif - - | threshold - -abs 0.7 $OUTDIR/${DWIFILENAME}_sf.mif

## Estimates the response function
estimate_response $OUTDIR/${DWIFILENAME}_dwi.mif $OUTDIR/${DWIFILENAME}_sf.mif -lmax 6 -grad $OUTDIR/$DWIFILENAME.b $OUTDIR/${DWIFILENAME}_response.txt

## Perform CSD in each white matter voxel
for i_lmax in 2 4 6 8 10 12; do
    csdeconv $OUTDIR/${DWIFILENAME}_dwi.mif -grad $OUTDIR/$DWIFILENAME.b $OUTDIR/${DWIFILENAME}_response.txt -lmax $i_lmax -mask $OUTDIR/${DWIFILENAME}_brainmask.mif $OUTDIR/${DWIFILENAME}_lmax${i_lmax}.mif
## echo DONE Lmax=$i_lmax 
done 

##
echo DONE performing preprocessing of data before starting tracking...
##

##
echo START tracking...
##
##echo tracking Deterministic Tensor based
for i_track in 01; do
streamtrack DT_STREAM $OUTDIR/${DWIFILENAME}_dwi.mif $OUTDIR/${DWIFILENAME}_wm_tensor-NUM${i_track}-$NUMFIBERS.tck -seed $OUTDIR/${DWIFILENAME}_wm.mif -mask $OUTDIR/${DWIFILENAME}_wm.mif -grad $OUTDIR/${DWIFILENAME}.b -number $NUMFIBERS -maxnum $MAXNUMFIBERSATTEMPTED
done

for i_track in 01; do
## loop over tracking and lmax
for i_tracktype in SD_STREAM SD_PROB; do
##
##echo Tracking $i_tracktype Deterministic=1 Probabilistic=2 CSD-based
##
    for i_lmax in 2 4 6 8 10 12; do
	##echo Tracking CSD-based Lmax=$i_lmax
	streamtrack $i_tracktype $OUTDIR/${DWIFILENAME}_lmax${i_lmax}.mif   $OUTDIR/${DWIFILENAME}_csd_lmax${i_lmax}_wm_${i_tracktype}-NUM${i_track}-$NUMFIBERS.tck -seed $OUTDIR/${DWIFILENAME}_wm.mif  -mask $OUTDIR/${DWIFILENAME}_wm.mif  -grad $OUTDIR/$DWIFILENAME.b -number $NUMFIBERS -maxnum $MAXNUMFIBERSATTEMPTED
    done
##
##echo DONE Tracking $i_tracktype Deterministic=1 Probabilistic=2 CSD-based
##
done
done
##
echo DONE tracking. Exiting Ensemble Tracking Candidate Fascicle Generation Script
##

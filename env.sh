# ~/redux/daily/run/scripts/night/20210322/poststdstar-20210322-00081536-a0123456789-40965715.log
# 
# grep RUNNING /global/homes/j/jguy/redux/daily/run/scripts/night/20210322/poststdstar-20210322-00081536-a0123456789-40965715.log | grep b0
# 
export SPECPROD=daily
export DESI_SPECTRO_CALIB=/global/cfs/cdirs/desi/spectro/desi_spectro_calib/trunk/  
export DESIMODEL=/global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/code/desimodel/master
export DESI_LOGLEVEL=DEBUG

#export DESI_SPECTRO_DATA=$CSCRATCH/trash/etc/desi/spectro/data

# desi_compute_fluxcalibration --infile /global/cfs/cdirs/desi/spectro/redux/daily/exposures/20210322/00081536/frame-b0-00081536.fits --sky /global/cfs/cdirs/desi/spectro/redux/daily/exposures/20210322/00081536/sky-b0-00081536.fits --fiberflat /global/cfs/cdirs/desi/spectro/desi_spectro_calib/trunk/spec/sm4/fiberflatnight-b0-20201214.fits --models /global/cfs/cdirs/desi/spectro/redux/daily/exposures/20210322/00081536/stdstars-0-00081536.fits --delta-color-cut 0.1 --outfile $CSCRATCH/trash/denali/fluxcalib-b0-00081536.fits

# desi_process_exposure --infile /global/cfs/cdirs/desi/spectro/redux/daily/exposures/20210322/00081536/frame-b0-00081536.fits --fiberflat /global/cfs/cdirs/desi/spectro/desi_spectro_calib/trunk/spec/sm4/fiberflatnight-b0-20201214.fits --sky /global/cfs/cdirs/desi/spectro/redux/daily/exposures/20210322/00081536/sky-b0-00081536.fits --calib $CSCRATCH/trash/denali/fluxcalib-b0-00081536.fits --cosmics-nsig 6 --outfile $CSCRATCH/trash/denali/cframe-b0-00081536.fits

# srun -N 16 -n 16 -c 16 desi_tsnr_afterburner --prod $SPECPROD --outfile $CSCRATCH/trash/denali/tsnr-exposures-$SPECPROD.fits --recompute --nproc 16 --aux /global/cfs/cdirs/desi/survey/observations/SV1/sv1-tiles.fits
# --nights 20200224,20200227,20200303,20200306,20201214,20201217,20201220,20201223,20210103,20210107,20210110,20210115,20210131,20210203,20210206,20210214,20210217,20210220,20210223

# srun -N 1 -n 1 -c 16 desi_tsnr_afterburner --prod $SPECPROD --outfile $CSCRATCH/trash/denali/tsnr-exposures-$SPECPROD.fits --recompute --nproc 16 --aux /global/cfs/cdirs/desi/survey/observations/SV1/sv1-tiles.fits --nights 20210224 --expids 77930

# desi_tsnr_afterburner --prod $SPECPROD --outfile $CSCRATCH/trash/tsnr-fiberloss/tsnr-exposures-$SPECPROD.fits --recompute --nproc 16 # --nights 20210321

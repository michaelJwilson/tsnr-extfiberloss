import os,sys
import glob
import yaml
import itertools
import argparse
import astropy.io.fits as fits
import fitsio
import numpy as np
import pylab as pl
import multiprocessing

from   astropy.table import Table,vstack
from   pkg_resources import resource_filename
from   desispec.io import read_sky
from   desispec.io import read_fiberflat
from   pathlib import Path
from   desispec.io.meta import findfile, specprod_root
from   desispec.calibfinder import CalibFinder
from   desispec.io import read_frame
from   desispec.io import read_fibermap
from   desispec.io.fluxcalibration import read_flux_calibration
from   desiutil.log import get_logger
from   desispec.tsnr import calc_tsnr2,tsnr2_to_efftime
from   astropy.table import Table, vstack
from   desiutil.depend import getdep
from   desispec.tilecompleteness import compute_tile_completeness_table,merge_tile_completeness_table
from   desispec.skymag import compute_skymag
from   desispec.efftime import compute_efftime

def compute_tsnr_values(cframe_filename,cframe_hdulist,night,expid,camera,specprod_dir, alpha_only=False) :
    """                                                                                                                                                                                                                          
    Computes TSNR values                                                                                                                                                                                                          
    Args:                                                                                                                                                                                                                        
       cframe_filename: str, cframe file path                                                                                                                                                                                    
       cframe_hdulist: astropy.fits.HDUlist object                                                                                                                                                                               
       night: int                                                                                                                                                                                                                
       expid: int                                                                                                                                                                                                                
       camera: str                                                                                                                                                                                                               
       specprod_dir: str, production directory                                                                                                                                                                                   
       alpha_only: bool, set to True to only recompute alpha                                                                                                                                                                                                                                                                                                                                                                                                     
    Returns: astropy.table.Table obkect with TSNR values                                                                                                                                                                        
    """

    calib  = findfile('fluxcalib', night=night, expid=expid,
                      camera=camera, specprod_dir=specprod_dir)
    flat = cframe_hdulist[0].header['FIBERFLT']
    if 'SPECPROD' in flat:
        flat = flat.replace('SPECPROD', specprod_dir)
    elif 'SPCALIB' in flat:
        hdr  = fitsio.read_header(cframe_filename)
        flat = flat.replace('SPCALIB', getdep(hdr, 'DESI_SPECTRO_CALIB'))
    else:
        raise ValueError('Failed on flat retrieval for {}.'.format(hdr))

    iin = cframe_filename.replace('cframe', 'frame')
    sky = cframe_filename.replace('cframe', 'sky')
    psf = cframe_filename.replace('cframe', 'psf')

    frame=read_frame(iin, skip_resolution=True)
    fiberflat=read_fiberflat(flat)
    fluxcalib=read_flux_calibration(calib)
    skymodel=read_sky(sky)

    results, alpha = calc_tsnr2(frame, fiberflat=fiberflat,
                                skymodel=skymodel, fluxcalib=fluxcalib, alpha_only=alpha_only)

    table=Table()
    for k in results:
        table[k] = results[k].astype(np.float32)
    table["TSNR2_ALPHA_"+camera[0].upper()] = np.repeat(alpha,len(frame.flux))

    return table


if __name__ == '__main__':
    night           = '20210319'
    expid           = '00081093'
    camera          = 'b1'

    specprod_dir    = '/project/projectdirs/desi/spectro/redux/denali/'
    
    cframe_filename = '{}/exposures/{}/{}/cframe-{}-{}.fits'.format(specprod_dir, night, expid, camera, expid)
    cframe_hdulist  = fits.open(cframe_filename) 

    print('\n\nDone.\n\n')

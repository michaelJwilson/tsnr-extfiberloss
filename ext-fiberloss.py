import os,sys
import glob
import yaml
import itertools
import argparse
import astropy.io.fits as fits
import fitsio
import pickle
import itertools
import numpy as np
import pylab as pl
import multiprocessing
import astropy.io.fits as fits

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

np.random.seed(314)

def compute_tsnr_values(cframe_filename,night,expid,camera,specprod_dir, alpha_only=False, model_extfiberloss=False, components=True, no_offsets=False) :
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

    cframe_hdulist = fits.open(cframe_filename)
    
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

    if components:
        results, alpha, seeing_fwhm, components = calc_tsnr2(frame, fiberflat=fiberflat,
                                                             skymodel=skymodel, fluxcalib=fluxcalib, alpha_only=alpha_only, model_extfiberloss=model_extfiberloss, components=True, no_offsets=no_offsets)
    else:
        results, alpha, seeing_fwhm = calc_tsnr2(frame, fiberflat=fiberflat,
                                                 skymodel=skymodel, fluxcalib=fluxcalib, alpha_only=alpha_only, model_extfiberloss=model_extfiberloss, no_offsets=no_offsets)
    
    table=Table()
    for k in results:
        table[k] = results[k].astype(np.float32)

    for k in components:
        table[k] = components[k].astype(np.float32)
        
    table["TSNR2_ALPHA_"+camera[0].upper()] = np.repeat(alpha,len(frame.flux))
    table["TSNR2_SEEING_"+camera[0].upper()] = np.repeat(seeing_fwhm,len(frame.flux)) 
    
    return table

def wrapper(args):
    return compute_tsnr_values(**args)

if __name__ == '__main__':
    arms            =  ['b','r','z']
    arms            =  ['r']

    petals          =  np.arange(0,10,1).astype(str)

    cameras         =  np.array([a+b for a in arms for b in petals])

    specprod_dir    = '/project/projectdirs/desi/spectro/redux/daily/'

    model_extfiberloss = True
    no_offsets         = True
    
    # exposures     = Table.read('/project/projectdirs/desi/survey/observations/SV1/sv1-exposures.fits')
    # exposures       = Table.read('/project/projectdirs/desi/spectro/redux/daily/tsnr-exposures.fits', 'TSNR2_EXPID')
    
    # exposures     = exposures[exposures['NIGHT'] == '20210324']
    # exposures     = Table(exposures[exposures['EXPID'] == 86986])

    # print('20210510', np.count_nonzero(exposures['NIGHT'] == '20210510'))
    
    # exposures     = Table(exposures[exposures['NIGHT'] == '20210510'])
    
    # exposures     = exposures[exposures['SURVEY'] == 'sv3']
    # exposures     = exposures[exposures['FAFLAVOR'] == 'sv3bright']
    # exposures     = exposures[exposures['TARGETS'] == 'BGS+MWS']
    # exposures     = exposures[exposures['GFA_FWHM_ASEC'] >= 0.0]
    # exposures.pprint()
    
    '''
    for x in exposures.dtype.names:
        print(x)
    ''' 
    # print(len(exposures))
    
    # rows            = np.random.randint(0, high=len(exposures), size=1, dtype=int)
    
    args = []
    auxs = []

    '''
    for row in rows:
        night            = exposures['NIGHT'][row]
        expid            = exposures['EXPID'][row]
        
        efftime_dark     = exposures['EFFTIME_DARK_GFA'][row]
        efftime_bright   = exposures['EFFTIME_BRIGHT_GFA'][row]

        seeing_etc       = exposures['SEEING_ETC'][row]
        efftime_etc      = exposures['EFFTIME_ETC'][row]
        
        camera           = np.random.choice(cameras, replace=True, size=1)[0]
        exp_seeing_fwhm  = exposures['SEEING_GFA'][row]        
        cframe_filename  = '{}/exposures/{}/{:08d}/cframe-{}-{:08d}.fits'.format(specprod_dir, night, expid, camera, expid)

        print(expid, night, camera, exp_seeing_fwhm, cframe_filename)

        if os.path.exists(cframe_filename):
            auxs.append({'row': row, 'efftime_dark': efftime_dark, 'efftime_bright': efftime_bright, 'seeing_fwhm': exp_seeing_fwhm, 'expid': expid, 'seeing_etc': seeing_etc, 'efftime_etc': efftime_etc})
            args.append({'cframe_filename':cframe_filename,'night':night,'expid':expid,'camera':camera,\
                         'specprod_dir':specprod_dir,'alpha_only':False})

        else:
            print('Failed to find: {}'.format(cframe_filename))
            continue
    '''

    files = glob.glob('/project/projectdirs/desi/spectro/redux/daily/exposures/20210510/*/cframe-r?-*.fits')

    for row, x in enumerate(files):
        parts = x.split('/')

        night = np.int(parts[-3])
        expid = np.int(parts[-2])

        camera = x.split('-')[-2]

        if camera != 'r0':
            print('skippind {} of {}'.format(camera, expid))
            
            continue

        cframe_filename  = '{}/exposures/{}/{:08d}/cframe-{}-{:08d}.fits'.format(specprod_dir, night, expid, camera, expid)

        meta = fits.open(cframe_filename)[0].header
        
        auxs.append({'row': row, 'night': night, 'expid': expid, 'camera': camera, 'fibassgb': meta['FIBASSGN'], 'program': meta['PROGRAM'], 'tileid': meta['TILEID']})
        args.append({'cframe_filename':cframe_filename,'night':night,'expid':expid,'camera':camera, 'specprod_dir':specprod_dir,'alpha_only':False , 'components': True,\
                     'model_extfiberloss': model_extfiberloss, 'no_offsets': no_offsets})     
        
        print(night, expid, camera, cframe_filename)
        
    with multiprocessing.Pool(8) as pool:
        results = pool.map(wrapper, args)
    
    for aux, arg, result in zip(auxs, args, results):
        aux.update(arg)

        aux.update({'tsnr_seeing': result['TSNR2_SEEING_R'][0]})
        aux.update({'tsnr_alpha':  result['TSNR2_ALPHA_R'][0]})

        for band in ['R']:
            for tt in ['ELG', 'BGS', 'QSO']:
                aux.update({'tsnr2_{}_{}'.format(tt.lower(), band.lower()): np.median(result['TSNR2_{}_{}'.format(tt, band)])})
                
                aux.update({'tsnr2_{}_{}_signal'.format(tt.lower(), band.lower()): np.median(result['TSNR2_{}_{}_SIGNAL'.format(tt, band)])})
                aux.update({'tsnr2_{}_{}_background'.format(tt.lower(), band.lower()): np.median(result['TSNR2_{}_{}_BACKGROUND'.format(tt, band)])})
                
    pickle.dump(auxs, open('/global/cscratch1/sd/mjwilson/trash/test_20210510_extfiberloss_{:d}_nooffsets_{:d}.pickle'.format(np.int(model_extfiberloss), np.int(no_offsets)), 'wb'))
        
    # pl.plot(aux['efftime_bright'], np.mean(result['TSNR2_BGS_R']), marker='.', lw=0.0, c='k', markersize=5)
    # pl.show()

    print('\n\nDone.\n\n')

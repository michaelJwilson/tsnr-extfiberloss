import os
import glob
import fitsio
import scipy
import numpy as np
import pylab as pl

import astropy.io.fits as fits

from   astropy import units, constants
from   desispec.io import read_sky, read_average_flux_calibration, findfile
from   desispec.io.fluxcalibration import read_flux_calibration
from   desispec.skymag import compute_skymag
from   speclite import filters
from   astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel


night=np.int('20210510')  # expid=np.int('87942')

fiber=0

specprod_dir = '/project/projectdirs/desi/spectro/redux/daily/'

expids = glob.glob('/project/projectdirs/desi/spectro/redux/daily//exposures/20210510/*')
expids = np.array([x.split('/')[-1] for x in expids]).astype(np.int)

wmin, wmax, wdelta = 3600, 9824, 0.8
fullwave = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
cslice = {"b": slice(0, 2751), "r": slice(2700, 5026), "z": slice(4900, 7781)}

# "decam2014-g", "decam2014-r", "decam2014-z"                                                                                                                                                                                              
filts = filters.load_filters("decam2014-g", "decam2014-r", "decam2014-z")

# 0.8 x npix.
gauss_kernel = Gaussian1DKernel(25)

for expid in expids:
    csky = np.zeros(fullwave.shape)
    ok   = True

    counts={}
    
    for band in ['b','r','z']:
        cframe='/project/projectdirs/desi/spectro/redux/daily//exposures/{}/{:08d}/cframe-{}0-{:08d}.fits'.format(night, expid, band, expid)
        sky=cframe.replace('cframe', 'sky')
    
        filename = findfile("sky",night=night,expid=expid,camera='{}0'.format(band),specprod_dir=specprod_dir)

        # print(sky)
        # print(filename)

        if not os.path.exists(filename):
            ok=False
            
            break
        
        header=fitsio.read_header(sky)
        exptime=header["EXPTIME"]

        cal_filename="{}/spec/fluxcalib/fluxcalibnight-{}-20201216.fits".format(os.environ["DESI_SPECTRO_CALIB"],band)
        
        fiber_acceptance_for_point_sources = 0.60 # see DESI-6043
        mean_fiber_diameter_arcsec = 1.52 # see DESI-6043
        fiber_area_arcsec = np.pi*(mean_fiber_diameter_arcsec/2)**2
    
        acal = read_average_flux_calibration(cal_filename)

        # calib=findfile('fluxcalib', night=night, expid=expid, camera='{}0'.format(band), specprod_dir=specprod_dir)
        # acal=read_flux_calibration(calib)
        
        # [e/A]
        sky=read_sky(sky)

        csky[cslice[band]] = sky.flux[fiber] / exptime / acal.value() * fiber_acceptance_for_point_sources / fiber_area_arcsec * 1e-17 # [ergs/s/cm2/A/arcsec2]

        # smooth_sky = convolve(sky.flux[fiber], gauss_kernel)
        smooth_sky   = csky[cslice[band]]  # scipy.signal.medfilt(csky[cslice[band]], 25)
        
        counts[band] = np.sum(sky.wave * smooth_sky)

        # pl.plot(sky.wave, sky.flux[fiber])
        # pl.plot(sky.wave, smooth_sky)
        # pl.title('{} {}'.format(expid, '{}0'.format(band)))
        # pl.show()
        
    if not ok:
        continue
        
    sky_pad, fullwave_pad = csky.copy(), fullwave.copy()
    
    for i in range(len(filts)):
        sky_pad, fullwave_pad = filts[i].pad_spectrum(sky_pad, fullwave_pad, method="zero")
        
    mags = filts.get_ab_magnitudes(sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom),fullwave_pad * units.angstrom).as_array()[0]
    smag = compute_skymag(night, expid, specprod_dir=specprod_dir)
    
    print('---- {:d} ----'.format(expid))
    print(smag)
    print(mags)
    print(counts)

    pl.plot(10.**(-mags[1] / 2.5), counts['r'], marker='.')

pl.show()

import pickle
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from   astropy.table import Table
from   desimodel.fastfiberacceptance import FastFiberAcceptance
from   desimodel.io import load_platescale

fa = FastFiberAcceptance()
ps = load_platescale()
isotropic_platescale = np.interp(np.array([0.0]),ps['radius']**2,np.sqrt(ps['radial_platescale']*ps['az_platescale'])) # um/arcsec 

def surveyspeed_fiberfrac(exposure_seeing_fwhm, ttype='bgs'):
    if ttype == 'bgs':
        return np.exp(0.0341 * np.log(exposure_seeing_fwhm)**3 -0.3611 * np.log(exposure_seeing_fwhm)**2 -0.7175 * np.log(exposure_seeing_fwhm) -1.5643)

    else:
        return np.exp(0.0989 * np.log(exposure_seeing_fwhm)**3 -0.5588 * np.log(exposure_seeing_fwhm)**2 -0.9708 * np.log(exposure_seeing_fwhm) -0.4473)
    
def fa_fiberfrac(exposure_seeing_fwhm, point=False):
    # we could include here a wavelength dependence on seeing, or non-Gaussian.                                                                                                                                                        
    sigmas_um = (exposure_seeing_fwhm/2.35) * isotropic_platescale # um   
    offsets_um = np.zeros_like(sigmas_um)

    if point:
        return  fa.value("POINT",sigmas_um)

    else:
        return  fa.value("BULGE",sigmas_um,1.50 * np.ones_like(sigmas_um))


if __name__ == '__main__':
    fwhms = np.arange(0.0, 10., 0.01)    

    fname = '/global/cscratch1/sd/mjwilson/trash/aux_psf.pickle'
    # fname = '/global/cscratch1/sd/mjwilson/trash/aux.pickle'                                                                                                                                                                      
    # fname = '/global/cscratch1/sd/mjwilson/trash/aux_poisson.pickle'                                                                                                                                                          
    # fname = '/global/cscratch1/sd/mjwilson/trash/aux_nooffsets.pickle'                                                                                                                                                                     
    auxs  = pickle.load( open( fname , "rb" ) )

    to_solve = []

    for aux in auxs:        
        to_solve.append([aux['efftime_bright'], aux['tsnr2_bgs_r'], aux['seeing_fwhm']])

    to_solve = np.array(to_solve)
    
    fig, axes = plt.subplots(1,3,figsize=(15,5))
    
    axes[0].hist(to_solve[:,2], bins=np.arange(0.0, 4.0, 0.05), histtype='step')

    #
    axes[1].plot(fwhms, surveyspeed_fiberfrac(fwhms), label='bgs survey speed')
    axes[1].plot(fwhms, fa_fiberfrac(fwhms), label='fa - zero offset')

    axes[1].set_xlabel('SEEING FWHM')
    axes[1].set_ylabel('FIBERFRAC')
    
    axes[1].legend(frameon=False)

    #
    axes[2].plot(fwhms, surveyspeed_fiberfrac(fwhms, ttype='backup'), label='bgs survey speed')
    axes[2].plot(fwhms, fa_fiberfrac(fwhms, point=True), label='fa - zero offset')

    axes[2].set_xlabel('SEEING FWHM')
    axes[2].set_ylabel('FIBERFRAC')

    axes[2].legend(frameon=False)
    
    # pl.show()

    pl.clf()

    exposures = Table.read('/project/projectdirs/desi/survey/observations/SV1/sv1-exposures.fits')
    exposures = exposures[exposures['TARGETS'] == 'BGS+MWS']
    exposures = exposures[exposures['GFA_FWHM_ASEC'] >= 0.0]
    exposures = exposures[exposures['GFA_FRACFLUX_NOMINAL_BGS'] > 0.0]
    exposures.sort('GFA_FWHM_ASEC')
    exposures.pprint()

    pl.plot(exposures['GFA_FWHM_ASEC'], exposures['GFA_FIBER_FRACFLUX_BGS'])
    pl.plot(exposures['GFA_FWHM_ASEC'], surveyspeed_fiberfrac(exposures['GFA_FWHM_ASEC'], ttype='bgs'), c='k')
    pl.xlabel('GFA_FWHM_ASEC')
    pl.ylabel('GFA_FIBER_FRACFLUX_BGS')
    pl.show()
    
    print('\nDone.\n')


    

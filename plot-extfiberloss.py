import json
import pickle
import scipy
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

from   scipy import stats
from   astropy.table import Table, vstack
from   desispec.io.meta import findfile

# fname = '/global/cscratch1/sd/mjwilson/trash/aux_psf.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux_poisson.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux_nooffsets.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux_nooffsets_nopoisson.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets_nopoisson.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets_nopoisson_noextended.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets_noextended.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_noextended.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent.pickle'
fname   = '/global/cscratch1/sd/mjwilson/trash/test_20210510.pickle'

auxs    = pickle.load( open( fname , "rb" ) ) 

to_solve = []

tables = []

camera = 'r0'

for aux in auxs:
    # print(aux['seeing_fwhm'])
    # to_solve.append([aux['efftime_etc'], aux['tsnr2_bgs_r'], aux['seeing_etc']])

    if aux['program'] == 'BRIGHT':
        aux['tsnr2_spec_r'] = aux['tsnr2_bgs_r']
    else:
        aux['tsnr2_spec_r'] = aux['tsnr2_elg_r']
        
    if aux['camera'] != camera:
        continue
    
    etcpath=findfile('etc', night=aux['night'], expid=aux['expid'])

    with open(etcpath) as f:
        etcdata = json.load(f)

    aux.update(etcdata['expinfo'])
    
    to_solve.append([aux['tsnr2_spec_r'], etcdata['expinfo']['efftime'], aux['tsnr2_bgs_r'], aux['tsnr2_elg_r']])  
    
    L1  = list(aux.keys())
    L2  = list(aux.values()) 
    
    aux = {k:[v] for k,v in zip(L1,L2)}
    
    tables.append(Table(aux)) 
#
#  'airmass': 1.181809, 'atm_extinction': 0.981091, 'realtime': 452.418579, 'signal': 0.432851, 'background': 1.364067, 'transp_obs': 0.708089, 'ffrac_psf': 0.308933, 'ffrac_elg': 0.249342, 'ffrac_bgs':, 'thru_psf': 0.221507
# 
tables = vstack(tables)
tables = tables['expid', 'tileid', 'night', 'program', 'airmass', 'atm_extinction', 'realtime', 'signal', 'background', 'transp_obs',\
                'ffrac_psf', 'ffrac_elg', 'ffrac_bgs', 'thru_psf', 'req_efftime', 'efftime', 'tsnr2_spec_r', 'tsnr2_bgs_r', 'tsnr2_elg_r',\
                'tsnr_seeing', 'tsnr_alpha']

tables.sort('tsnr2_spec_r')
tables.pprint()

# del tables['nominal_seeing_fwhm']

tables.write(fname.replace('.pickle', '.fits'), format='fits', overwrite=True)

fig, axes = plt.subplots(1, 2, figsize=(10,5))

for i, (marker, program) in enumerate(zip(['^', '*'], ['DARK', 'BRIGHT'])):
    isin = tables['program'] == program

    idx = np.argsort(tables['tsnr2_spec_r'][isin])
    
    res = stats.linregress(tables['tsnr2_spec_r'][isin][idx], tables['efftime'][isin][idx])

    axes[i].plot(tables['tsnr2_spec_r'][isin], res.intercept + res.slope * tables['tsnr2_spec_r'][isin], 'k', lw=0.5)
    
    axes[i].scatter(tables['tsnr2_spec_r'][isin], tables['efftime'][isin], marker=marker, lw=0.0, s=14)
    axes[i].set_xlabel('TSNR2_SPEC_R')
    axes[i].set_ylabel('efftime_etc')
    axes[i].set_title('{} {} (pearson r: {:.3f})'.format(program, camera, stats.pearsonr(tables['tsnr2_spec_r'][isin], tables['efftime'][isin])[0]))

    # pl.colorbar(label='SEEING FWHM')
    
pl.show() 

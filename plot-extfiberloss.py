import json
import pickle
import scipy
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

from   scipy import stats
from   astropy.table import Table, vstack, join
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

# fname = '/global/cscratch1/sd/mjwilson/trash/test_20210510_extfiberloss_1.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/test_20210510_extfiberloss_1_nooffsets_1_unitalpha.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/test_20210510_extfiberloss_1_nooffsets_1_sky_ivar_1_unitalpha_0.pickle'
fname = '/global/cscratch1/sd/mjwilson/trash/test_20210510_extfiberloss_1_nooffsets_1_sky_ivar_0_unitalpha_0.pickle'

auxs    = pickle.load( open( fname , "rb" ) ) 

tables  = []
camera  = 'r0'

for aux in auxs:    
    if aux['program'] == 'BRIGHT':
        aux['tsnr2_spec_r'] = np.array(aux['tsnr2_bgs_r'], copy=True)
    else:
        aux['tsnr2_spec_r'] = np.array(aux['tsnr2_elg_r'], copy=True)
        
    if aux['camera'] != camera:
        continue

    try:
        etcpath=findfile('etc', night=aux['night'], expid=aux['expid'])

        with open(etcpath) as f:
            etcdata = json.load(f)
    except:
        print('failed to find {}'.format(etcpath))
        continue
            
    aux.update(etcdata['expinfo'])
        
    L1 = list(aux.keys())
    L2 = [aux[x] for x in L1] 

    aux = {k:[v] for k,v in zip(L1,L2)}
    
    tables.append(Table(aux))
    
# ['airmass': 1.181809, 'atm_extinction': 0.981091, 'realtime': 452.418579, 'signal': 0.432851, 'background': 1.364067, 'transp_obs': 0.708089, 'ffrac_psf': 0.308933, 'ffrac_elg': 0.249342, 'ffrac_bgs':, 'thru_psf': 0.22150].
tables = vstack(tables)
'''
for x in tables.dtype.names:
    print(x)
'''

tables = tables['expid', 'tileid', 'night', 'program', 'sbprof', 'airmass', 'atm_extinction', 'realtime', 'signal', 'background', 'transp_obs',\
                'ffrac_psf', 'ffrac_elg', 'ffrac_bgs', 'thru_psf', 'req_efftime', 'efftime', 'tsnr2_spec_r', 'tsnr2_bgs_r', 'tsnr2_elg_r',\
                'tsnr_seeing', 'tsnr_alpha', 'tsnr2_elg_r_signal', 'tsnr2_elg_r_background', 'tsnr2_qso_r_signal', 'tsnr2_qso_r_background',\
                'tsnr2_bgs_r_background', 'tsnr2_bgs_r_signal']

tables.sort('tsnr2_spec_r')

daily  = Table.read('/project/projectdirs/desi/spectro/redux/daily/tsnr-exposures.fits', 1)
daily  = Table(daily[np.isin(daily['NIGHT'], np.array([20210510, 20210511, 20210512]))])
daily  = daily['EXPID', 'SKY_MAG_AB_GFA', 'SKY_MAG_G_SPEC', 'SKY_MAG_R_SPEC', 'SKY_MAG_Z_SPEC', 'EFFTIME_ETC', 'EFFTIME_SPEC', 'EFFTIME_GFA', 'EFFTIME_DARK_GFA', 'EFFTIME_BRIGHT_GFA', 'EFFTIME_BACKUP_GFA']

keys   = list(daily.dtype.names) 

for x in keys:
    daily[x.lower()] = daily[x.upper()]
    
    del daily[x.upper()]

daily.pprint()
    
tables = join(tables, daily, join_type='left', keys='expid')
tables.pprint()

tables.write(fname.replace('.pickle', '.fits'), format='fits', overwrite=True)

print('Writing: {}'.format(fname.replace('.pickle', '.fits')))


fig, axes = plt.subplots(1, 2, figsize=(10,5))

for i, (marker, program) in enumerate(zip(['^', '*'], ['DARK', 'BRIGHT'])):
    for color, night in zip(['orange','dodgerblue'],['20210510','20210511', '20210512']):
        isin = (tables['program'] == program) & (tables['night'] == night)
    
        gradient=np.sum(tables['tsnr2_spec_r'][isin] * tables['efftime_gfa'][isin]) / np.sum(tables['tsnr2_spec_r'][isin]**2.)
    
        # axes[i].scatter(tables['tsnr2_spec_r'][isin], tables['efftime'][isin], marker=marker, lw=0.0, s=14)
        # axes[i].set_ylabel('efftime_etc / TSNR2_SPEC_R')
    
        axes[i].scatter(tables['tsnr2_spec_r'][isin], tables['efftime'][isin] / (gradient * tables['tsnr2_spec_r'][isin]), marker=marker, lw=0.0, s=14, label=night)

    axes[i].set_xlabel('TSNR2_SPEC_R')
    axes[i].set_ylabel('efftime_etc / (TSNR2_SPEC_R | EFFTIME_GFA)')
    axes[i].set_title('{} {} (pearson r: {:.3f})'.format(program, camera, stats.pearsonr(tables['tsnr2_spec_r'][isin], tables['efftime'][isin])[0]))

pl.legend(frameon=False)
pl.show() 
pl.clf()

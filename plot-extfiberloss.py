import pickle
import scipy
import pylab as pl
import numpy as np

from   scipy import stats
from   astropy.table import Table, vstack

# fname = '/global/cscratch1/sd/mjwilson/trash/aux_psf.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux_poisson.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux_nooffsets.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux_nooffsets_nopoisson.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets_nopoisson.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets_nopoisson_noextended.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets_noextended.pickle'
fname = '/global/cscratch1/sd/mjwilson/trash/onepercent_nooffsets.pickle'
# fname   = '/global/cscratch1/sd/mjwilson/trash/onepercent_noextended.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/onepercent.pickle'

auxs    = pickle.load( open( fname , "rb" ) ) 

to_solve = []

tables = []

for aux in auxs:
    # print(aux['seeing_fwhm'])
    to_solve.append([aux['efftime_etc'], aux['tsnr2_bgs_r'], aux['seeing_etc']])

    L1  = list(aux.keys())
    L2  = list(aux.values()) 
    
    aux = {k:[v] for k,v in zip(L1,L2)}
    
    tables.append(Table(aux))
    
to_solve       =  np.array(to_solve)
to_solve[:,1] /= np.median(to_solve[:,1])

tables = vstack(tables)
tables.pprint()

del tables['nominal_seeing_fwhm']

tables.write(fname.replace('.pickle', '.fits'), format='fits', overwrite=True)

pl.scatter(to_solve[:,0], to_solve[:,1], c=to_solve[:,2], marker='.', lw=0.0, s=3, vmin=0.8, vmax=2.)

pl.xlabel('efftime_etc')
pl.ylabel('TSNR2_BGS_R / <TSNR2_BGS_R>')
pl.title('{} (pearson r: {:.3f})'.format(fname.split('/')[-1], stats.pearsonr(to_solve[:,0], to_solve[:,1])[0]))
pl.colorbar(label='SEEING FWHM')
pl.show() 

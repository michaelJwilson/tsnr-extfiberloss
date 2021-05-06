import pickle
import scipy
import pylab as pl
import numpy as np

from scipy import stats

# fname = '/global/cscratch1/sd/mjwilson/trash/aux_psf.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux.pickle'
# fname = '/global/cscratch1/sd/mjwilson/trash/aux_poisson.pickle'
fname = '/global/cscratch1/sd/mjwilson/trash/aux_nooffsets.pickle'

auxs  = pickle.load( open( fname , "rb" ) ) 

to_solve = []

for aux in auxs:
    # print(aux['seeing_fwhm'])
    to_solve.append([aux['efftime_bright'], aux['tsnr2_bgs_r'], aux['seeing_fwhm']])

to_solve       =  np.array(to_solve)
to_solve[:,1] /= np.median(to_solve[:,1])

pl.scatter(to_solve[:,0], to_solve[:,1], c=to_solve[:,2], marker='.', lw=0.0, s=3, vmin=0.8, vmax=2.)

pl.xlabel('efftime_bright')
pl.ylabel('TSNR2_BGS_R / <TSNR2_BGS_R>')
pl.title('{} (pearson r: {:.3f})'.format(fname.split('/')[-1], stats.pearsonr(to_solve[:,0], to_solve[:,1])[0]))
pl.colorbar(label='SEEING FWHM')
pl.show() 

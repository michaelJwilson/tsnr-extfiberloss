import os
import glob
import numpy as np


files = glob.glob('/project/projectdirs/desi/spectro/redux/daily/exposures/20210510/*/frame-z?-*.fits')

specprod_dir = '/project/projectdirs/desi/spectro/redux/daily/'

for x in files:
    parts = x.split('/')

    night = np.int(parts[-3])
    expid = np.int(parts[-2])

    camera = x.split('-')[-2]
    
    cframe_filename  = '{}/exposures/{}/{:08d}/cframe-{}-{:08d}.fits'.format(specprod_dir, night, expid, camera, expid)
    
    print(night, expid, camera, specprod_dir)
    
    break

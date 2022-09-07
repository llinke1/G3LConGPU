from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.table import Table
import pandas as pd
import time
import pyccl as ccl
pd.options.display.max_columns = 500




filename="/vol/aibn238/data1/llinke/LowZ_IA/data/lowz12_v5_short.fits"
outputfilename_shapes="/vol/aibn238/data1/llinke/LowZ_IA/data/shapes"
outputfilename_lenses="/vol/aibn238/data1/llinke/LowZ_IA/data/lenses"

data = Table.read(filename, format='fits').to_pandas()

data = data[data['Z']>0]

# Singh+2015 - WMAP9 Cosmology
cosmo = ccl.Cosmology(Omega_c = 0.236, Omega_b = 0.046, h = 0.7, sigma8 = 0.817, n_s = 0.9646)

data['DA'] = ccl.angular_diameter_distance(cosmo, 1/(1+data['Z']))*cosmo['h'] # in Mpc/h


z_mask = (data['Z']>0.16) &(data['Z']<0.36)
e_mask = (data['e1']>-900) & (data['e2']>-900)  & ((data['e1']**2+data['e2']**2)**0.5<2)
mr_mask_lens = (data['Mr_ke0']>-900) & (data['Mr_ke0']<-22.05)
mr_mask_shapes = (data['Mr_ke0']>-900) & (data['Mr_ke0']<-22.02)

mu_mask = (data['Mu_ke0']>-900) & (data['Mu_ke0']<-15)
mz_mask = (data['Mz_ke0']>-900) & (data['Mz_ke0']<-15)

density = data[z_mask& e_mask & mr_mask_lens & mu_mask & mz_mask]

shapes = data[z_mask & e_mask & mr_mask_shapes & mu_mask & mz_mask]

north_shapes = shapes[(shapes['RA']>60) & (shapes['RA']<280)]
south_shapes = shapes[(shapes['RA']<60) | (shapes['RA']>280)]

north_density = density[(density['RA']>60) & (density['RA']<280)]
south_density = density[(density['RA']<60) | (density['RA']>280)]



RA_cen, DEC_cen = north_shapes['RA'].mean(), north_shapes['DEC'].mean()
print(RA_cen, DEC_cen)

north_shapes['X'] = np.deg2rad(north_shapes['RA']-RA_cen)*north_shapes['DA']/np.cos(np.deg2rad(north_shapes['DEC']))
north_shapes['Y'] = np.deg2rad(north_shapes['DEC']-DEC_cen)*north_shapes['DA']

north_density['X'] = np.deg2rad(north_density['RA']-RA_cen)*north_density['DA']/np.cos(np.deg2rad(north_density['DEC']))
north_density['Y'] = np.deg2rad(north_density['DEC']-DEC_cen)*north_density['DA']


RA_cen, DEC_cen = 0., south_shapes['DEC'].mean()
print(RA_cen, DEC_cen)

south_shapes['X'] = np.deg2rad(south_shapes['RA']-RA_cen)*south_shapes['DA']/np.cos(np.deg2rad(south_shapes['DEC']))
south_shapes['Y'] = np.deg2rad(south_shapes['DEC']-DEC_cen)*south_shapes['DA']

south_density['X'] = np.deg2rad(south_density['RA']-RA_cen)*south_density['DA']/np.cos(np.deg2rad(south_density['DEC']))
south_density['Y'] = np.deg2rad(south_density['DEC']-DEC_cen)*south_density['DA']


column_list=['X', 'Y', 'e1', 'e2', 'Z', 'WEIGHT_SEEING']


north_shapes[column_list].to_csv(outputfilename_shapes+"_n.dat",index=False, sep=' ')
south_shapes[column_list].to_csv(outputfilename_shapes+"_s.dat",index=False, sep=' ')

north_density[column_list].to_csv(outputfilename_lenses+"_n.dat",index=False, sep=' ')
south_density[column_list].to_csv(outputfilename_lenses+"_s.dat",index=False, sep=' ')

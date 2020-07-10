#!/usr/bin/env python
# coding: utf-8

# In[17]:


import os
from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import h5py

import math
import ephem
from datetime import datetime, time, timedelta
from math import pi, cos, sin
from scipy import optimize
from joblib import Parallel, delayed
from mpl_toolkits.basemap import Basemap
# import xesmf as xr

from diurnal_utils import *
# %run diurnal_utils.py


# In[18]:


start_date = '2000-06'
end_date =  '2010-06'
# end_date =  '2001-06'


# In[19]:


save_output_dir = '/export/data1/cchristo/diurnal_analysis_results/GPM/'
save_figs_dir = '/home/cchristo/proj_tapio/figs/diurnal_cycle_figs/gpm_2000_2010_precip/'


# In[57]:


# model_dir = '/export/data1/cchristo/gpm_data/gpmdata/'
model_dir = '/export/data1/cchristo/gpm_data/gpmdata_regridded_gfdl_cm4/'

years_to_include = range(2000, 2012)
years_to_include = [str(year_ii) for year_ii in years_to_include]

result = [y for x in os.walk(model_dir) for y in glob(os.path.join(x[0], '*.nc'))]


# filter paths based on year
result_filt =[]
for path_i in result: 
    if (os.path.dirname(path_i).split('/')[-4] in years_to_include):
        result_filt.append(path_i)

# get only certain years
# result = sorted(result)
result = sorted(result_filt)


# In[ ]:





# In[5]:


# %%time
ds = xr.open_mfdataset(result, combine='nested', concat_dim = 'time')
# ds = xr.open_dataset(result[0],  engine = 'h5netcdf')


# In[82]:


# ds['lon'] = ds['lon'].where(ds['lon'] > 0, ds['lon'] + 360)


# In[ ]:





# In[6]:


ds = ds.sel(time = slice(start_date, end_date))
ds = ds.sel(lat= slice(-60, 60))
# ds.load()
# print('Dataset Loaded')


# In[10]:


# ds['precipitationCal'].fillna(0)


# In[ ]:


# ds_mean = ds['precipitationCal'].mean(dim = 'time')
# ds_mean.load()

# ds_max = ds['precipitationCal'].mean(dim = 'time')
# ds_max.load()


# In[ ]:


# ds_mean.plot.imshow(vmin =0, vmax = 1, cmap = plt.get_cmap('gist_ncar'))
# ds_max.plot.imshow(vmin =0, vmax = 0.25 , cmap = plt.get_cmap('gist_ncar'))


# In[ ]:


# ds_mean.plot.imshow(vmin =0, vmax = 2, cmap = plt.get_cmap('gist_ncar'))


# In[ ]:


# ds.isel(time = 0)['precipitationCal'].plot.imshow(vmin =0, vmax = 10)


# In[ ]:


# ds['precipitationCal']


# In[8]:



out_ds, out_ds_means = diurnal_analysis(ds, 
                                        field_id = 'precipitationCal', 
                                        grid_time_resolution_hours = 0.5,
                                        time_resolution_hours = 0.5)


# In[9]:


out_ds.to_netcdf(save_output_dir + start_date + '_' + end_date + '_precip.nc')
out_ds_means.to_netcdf(save_output_dir + start_date + '_' + end_date + '_precip_diurnal_means.nc')


# In[ ]:





# In[ ]:


# make_four_panel(mu_season , 
# #                 cmap = plt.get_cmap('bwr'),
#                 cmap = plt.get_cmap('gist_ncar'),
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
#                 vmin = 0,
#                 vmax = 0.8,
#                 title = r'$\mu$',
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_means_pr.png')


# In[ ]:





# In[ ]:


# make_four_panel(sigma_season , 
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
#                 vmin = 0, vmax = 0.2, 
# #                 cmap = plt.get_cmap('bwr'),
#                 cmap = plt.get_cmap('gist_ncar'),
#                 title = r'$\sigma$',
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_stds_pr.png')


# In[ ]:


# make_four_panel(ampl_season , 
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
# #                 vmax = 0.00003, 
# #                 cmap = plt.get_cmap('bwr'),
#                 cmap = plt.get_cmap('gist_ncar'),
#                 title = r'$A$',
#                 vmin = 0, vmax = 0.2,
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_ampl_pr.png')


# In[ ]:


# make_four_panel(phase_season , 
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
# #                 vmin = 0, vmax = 24, 
#                 cmap = plt.get_cmap('twilight'),
#                 title = r'$\Phi$',
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_phase_pr.png')


# In[18]:


# out_ds = xr.Dataset()
# out_ds['mu_season'] = make_da_from_dict(mu_season, ds)
# out_ds['sigma_season'] = make_da_from_dict(sigma_season, ds)
# out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds)
# out_ds['phase_season'] = make_da_from_dict(phase_season, ds)


# In[22]:


# out_ds.to_netcdf(save_output_dir + 'gpm_2000_2010_precip.nc')


#!/usr/bin/env python
# coding: utf-8

# In[11]:


import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import ephem
from datetime import datetime, time, timedelta
from scipy import optimize
from mpl_toolkits.basemap import Basemap

import diurnal_config

from diurnal_utils import *
from fetch_model_helper import *
# %run fetch_model_helper.py


# # Notebook for exploring local CMIP6 data downloaded with `cmip6_downloader.py`

# In[12]:


field_name = 'clt'

cmip_identifier = 'CMIP6'
start_date = '1985-01'
end_date = '2006-01'

# cmip_identifier = 'CMIP5'
# start_date = '1985-01'
# end_date = '2006-01'



#TODO overwrite_existing_files = False


# In[13]:


# get all available model names
rel_cmip6_path= '/export/data1/cchristo/CMIP6_clouds/clt_3hr_historical/'
rel_cmip5_path = '/export/data1/cchristo/CMIP5_clouds/clt_3hr_historical/'
unique_cmip6_models = get_unique_models(rel_cmip6_path)
# unique_cmip5_models = get_unique_models(rel_cmip5_path)


# In[18]:




if cmip_identifier == 'CMIP6':
    all_model_names = unique_cmip6_models
    cmip_rel_dir = rel_cmip6_path

# elif cmip_identifier == 'CMIP5':
#     all_model_names = unique_cmip5_models
#     cmip_rel_dir = rel_cmip5_path
    


# In[16]:


# unique_cmip6_models


# In[ ]:



for model_name in list(all_model_names):
    print('Started... ', model_name)
    save_output_dir = '/export/data1/cchristo/diurnal_analysis_results/' + cmip_identifier + '_clouds/clt/'+ model_name + '/' 
    
    save_output_path = save_output_dir + start_date + '_' + end_date + '_' + field_name + '.nc'
    save_output_path_means = save_output_dir + start_date + '_' + end_date + '_' + field_name + '_diurnal_means.nc'
     # make dirs if they don't already exist

    if not os.path.exists(save_output_dir):
        os.makedirs(save_output_dir)
    
    # if files already exist, skip
    if (not os.path.exists(save_output_path)) &         (not os.path.exists(save_output_path_means)):
   
        try:
            #### Load data
            if cmip_identifier == 'CMIP6':
                path_to_cmip_files =  get_path_to_desired_model_cmip6(cmip_rel_dir, 
                                      desired_model= model_name,
                                      desired_ensemble_member = ('r1i1p1f1', 'r1i1p2f1','r1i1p1f2'),
                                      desired_grid_types = ('gn', 'gr', 'gr1', 'gr2'))
            elif cmip_identifier == 'CMIP5':
                path_to_cmip_files = get_path_to_desired_model_cmip5(cmip_rel_dir, 
                                  desired_model=model_name,
                                  desired_ensemble_member = ('r1i1p1','r6i1p1','r2i1p1'))
            # subset lat/lon and time
            print('Opening data...')
            ds = xr.open_mfdataset(path_to_cmip_files, combine='by_coords')
            ds = ds.sel(time = slice(start_date, end_date))
            ds = ds.sel(lat= slice(-60, 60))

            # perform diurnal analysis 
            print('Performing diurnal analysis... ')
            ds_sub = ds[field_name].to_dataset()

            out_ds, out_ds_means = diurnal_analysis(ds_sub, 
                                                    field_id = field_name, 
                                                    grid_time_resolution_hours = 3,
                                                    time_resolution_hours = 1)
#             # add some metadata 
            out_ds.attrs['input_dataset_paths'] = path_to_cmip_files
            out_ds_means.attrs['input_dataset_paths'] = path_to_cmip_files
            
            # save results 
            print('Saving results... ')
            out_ds.to_netcdf(save_output_path)
            out_ds_means.to_netcdf(save_output_path_means)
            
        except Exception as e:
            print('Could not process ' + model_name)
            print(e)

print('DONE!')


# # Fit diurnal cloud cycle with sin fit

# In[17]:


# mu_season, sigma_season, ampl_season, phase_season = diurnal_analysis(ds, 
#                                                                       field_id = 'clt', 
#                                                                       grid_time_resolution_hours=3,
#                                                                       time_resolution_hours = 1)


# In[30]:





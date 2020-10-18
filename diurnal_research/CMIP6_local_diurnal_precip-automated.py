#!/usr/bin/env python
# coding: utf-8

# In[14]:


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

# In[15]:


out_folder_identifier = 'bin_pt2_hr_round_gpm_overlap' #label to append to output path
skip_if_folder_exists =  True # used for processing models in parallel



# cmip_identifier = 'CMIP6'
# start_date = '1985-01'
# end_date = '2006-01'


cmip_identifier = 'CMIP6'
start_date = '2000-01'
end_date = '2015-01'

# cmip_identifier = 'CMIP5'
# start_date = '1985-01'
# end_date = '2006-01'

#TODO overwrite_existing_files = False


# In[16]:


# get all available model names
rel_cmip6_path= '/export/data1/cchristo/CMIP6_precip/pr_3hr_historical/'
rel_cmip5_path = '/export/data1/cchristo/CMIP5_precip/pr_3hr_historical/'
unique_cmip6_models = get_unique_models(rel_cmip6_path)
unique_cmip5_models = get_unique_models(rel_cmip5_path)


# In[35]:


# unique_cmip6_models
# unique_cmip6_models
# len(unique_cmip6_models)


# In[4]:




if cmip_identifier == 'CMIP6':
    all_model_names = unique_cmip6_models
    cmip_rel_dir = rel_cmip6_path

elif cmip_identifier == 'CMIP5':
    all_model_names = unique_cmip5_models
    cmip_rel_dir = rel_cmip5_path
    


# In[5]:


# get_path_to_desired_model_cmip6(cmip_rel_dir, 
#                                   desired_model= 'CanESM5',
#                                   desired_grid_types = ('gn', 'gr', 'gr1'))
# list(all_model_names)
# all_model_names = ['IPSL-CM6A-LR',]
# all_model_names


# In[7]:



for model_name in list(all_model_names):
    print('Started... ', model_name)

    save_output_dir = '/export/data1/cchristo/diurnal_analysis_results/' +         cmip_identifier + '_'+ out_folder_identifier + '/' + model_name + '/'
    
    save_output_path = save_output_dir + start_date + '_' + end_date + '_precip.nc'
    save_output_path_means = save_output_dir + start_date + '_' + end_date + '_precip_diurnal_means.nc'
    
    # skip folder if it's already been created (another script may be processing it)
    skip_model_iteration = False
    if skip_if_folder_exists & os.path.exists(save_output_dir):
        skip_model_iteration = True

    else:
         # make dirs if they don't already exist
        if not os.path.exists(save_output_dir):
            os.makedirs(save_output_dir)
    
    # if files already exist, skip
    if (not os.path.exists(save_output_path)) &         (not os.path.exists(save_output_path_means)) &         (not skip_model_iteration):
   
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
            ds_sub = ds['pr'].to_dataset()

            out_ds, out_ds_means = diurnal_analysis(ds_sub, 
                                                    field_id = 'pr', 
                                                    grid_time_resolution_hours = 3,
                                                    time_resolution_hours = 0.2)
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


# In[12]:


# ds.isel(time = 0)['time'].values.item()


# In[112]:


# %run fetch_model_helper.py
# path_to_cmip_files =  get_path_to_desired_model_cmip6(cmip_rel_dir, 
#                                       desired_model= 'EC-Earth3-Veg-LR', #model_name,
#                                       desired_ensemble_member = ('r1i1p1f1', 'r1i1p2f1','r1i1p1f2'),
#                                       desired_grid_types = ('gn', 'gr', 'gr1', 'gr2'))


# In[101]:


# path_to_cmip_files
# stacked = np.stack(list(average_cycle_season.values()))
# save_output_dir + start_date + '_' + end_date + '_precip.nc'


# In[102]:


# print(stacked.shape)
# # plt.plot(stacked[0,:,100,100])
# plt.plot(average_cycle_season['DJF'][:,50,50])


# In[32]:


# res.isel(season = 0, lat = 50, lon = 50).plot()


# In[ ]:





# In[35]:


# FLUX_TO_MM_HR = 60*60


# In[6]:


# model_list = os.listdir(path_to_cmip_dirs)
# print(model_list)
# path_to_cmip_files = path_to_cmip_dirs + 'GFDL-CM4/'
# file_list = os.listdir(path_to_cmip_files)
# for file in file_list: print(file)


# In[78]:



# ds = xr.open_mfdataset(path_to_cmip_files, combine='by_coords')

# # ds = ds.sel(time=slice('1986','2005'))
# # ds = ds.sel(time = slice('2000-06', '2001-05'))
# ds = ds.sel(time = slice('1999-01', '2015-01'))
# ds = ds.sel(lat= slice(-60, 60))


# In[77]:


# ds_sub = ds['pr'].to_dataset()
# mu_season, sigma_season, ampl_season, phase_season = diurnal_analysis(ds_sub, 
#                                                                       field_id = 'pr', 
#                                                                       grid_time_resolution_hours = 3,
#                                                                       time_resolution_hours = 1)


# In[76]:


# mu_mm_hr = {key:FLUX_TO_MM_HR*val for key, val in mu_season.items()}
# make_four_panel(mu_mm_hr, 
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
# #                 cmap = plt.get_cmap('bwr'),
#                 cmap = plt.get_cmap('gist_ncar'),
#                 vmin = 0,
#                 vmax = 0.8,
#                 title = r'$\mu$',
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_means_pr' + start_date + '_' + end_date +'.png')


# In[ ]:





# In[75]:


# sigma_mm_hr = {key:FLUX_TO_MM_HR*val for key, val in sigma_season.items()}


# make_four_panel(sigma_mm_hr , 
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
# #                 vmax = 0.1, 
#                 vmin = 0, vmax = 0.2, 
# #                 cmap = plt.get_cmap('bwr'),
#                 cmap = plt.get_cmap('gist_ncar'),
#                 title = r'$\sigma$',
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_stds_pr.png')


# In[74]:


# ampl_mm_hr = {key:FLUX_TO_MM_HR*val for key, val in ampl_season.items()}



# make_four_panel(ampl_mm_hr, 
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
# #                 vmax = 0.000015, 
# #                 cmap = plt.get_cmap('bwr'),
#                 vmin = 0, vmax = 0.2, 
#                 cmap = plt.get_cmap('gist_ncar'),
#                 title = r'$A$',
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_ampl_pr.png')


# In[73]:


# make_four_panel(phase_season , 
#                 lats = ds['lat'].values, 
#                 lons = ds['lon'].values,
#                 vmin = 0, vmax = 24, 
#                 cmap = plt.get_cmap('twilight'),
#                 title = r'$\Phi$',
# #                 axis = plt.axis([220, 300, 10, 50]), 
#                 save_fig_path= save_figs_dir + 'GFDL_CM4_phase_pr.png')


# In[71]:


# out_ds = xr.Dataset()
# out_ds['mu_season'] = make_da_from_dict(mu_season, ds)
# out_ds['sigma_season'] = make_da_from_dict(sigma_season,ds)
# out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds)
# out_ds['phase_season'] = make_da_from_dict(phase_season,ds)
# out_ds.to_netcdf(save_output_dir + 'gfdl_cm4_2000_2010_precip.nc')


# In[72]:


# out_ds.to_netcdf(save_output_dir + 'gfdl_cm4_2000_2010_precip.nc')


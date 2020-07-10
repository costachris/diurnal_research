#!/usr/bin/env python
# coding: utf-8

# ### Regridding
# Tools for regridding GPM data 
# from native satellite grid to model grid.

# In[38]:


import os
from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import xesmf as xe


# In[33]:


regrid_weights_dir = '/export/data1/cchristo/regridding_weights/'


# ## Load model grid

# In[12]:


path_to_cmip_dirs = '/export/data1/cchristo/cmip6_clouds/download_from_source/'
path_to_cmip_files = path_to_cmip_dirs + 'GFDL-CM4/'
file_list = os.listdir(path_to_cmip_files)

ds_model_grid = xr.open_dataset(path_to_cmip_files + file_list[2])


# In[29]:


grid_in = xr.Dataset({'lat': ds_model_grid['lat'].values,
                       'lon': ds_model_grid['lon'].values})


# ## Load Satellite Grid

# In[19]:


model_dir = '/export/data1/cchristo/gpm_data/gpmdata/'

result = [y for x in os.walk(model_dir) for y in glob(os.path.join(x[0], '*.nc'))]
result = sorted(result)

ds_satellite_grid = xr.open_dataset(result[10])


# In[30]:


grid_out = xr.Dataset({'lat': ds_satellite_grid['lat'].values,
                       'lon': ds_satellite_grid['lon'].values})


# In[32]:


grid_in


# In[31]:


grid_out


# In[39]:


regridder = xe.Regridder(grid_in, grid_out, 
                         method = 'bilinear', 
                         reuse_weights = False, 
                         filename = regrid_weights_dir + 'gpm_to_gfdl_cm4.nc')


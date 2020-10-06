#!/usr/bin/env python
# coding: utf-8

# In[16]:


import os
from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import xesmf as xe
import time
import gc
import grid_definitions

import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


# In[18]:


regrid_weights_dir = '/export/data1/cchristo/regridding_weights/'
regridded_output_dir = '/export/data1/cchristo/gpm_data/gpmdata_regridded_gfdl_cm4/'
# regridded_output_dir = '/export/data1/cchristo/gpm_data/gpmdata_regridded_grid1/'


# In[19]:


path_to_cmip_dirs = '/export/data1/cchristo/cmip6_clouds/download_from_source/'
path_to_cmip_files = path_to_cmip_dirs + 'GFDL-CM4/'
file_list = os.listdir(path_to_cmip_files)

ds_model_grid = xr.open_dataset(path_to_cmip_files + file_list[-1])


# In[21]:


grid_out = xr.Dataset({'lat': ds_model_grid['lat'].values,
                       'lon': ds_model_grid['lon'].values})


# In[27]:


# np.diff(grid_out.lon)


# In[53]:


# # grid_out.lon
# out_grid_name = 'grid1'
# grid_out = grid_definitions.grids[out_grid_name]


# In[28]:


model_dir = '/export/data1/cchristo/gpm_data/gpmdata/'

result = [y for x in os.walk(model_dir) for y in glob(os.path.join(x[0], '*.nc'))]
result = sorted(result)

ds_satellite_grid = xr.open_dataset(result[10])



# In[29]:


grid_in = xr.Dataset({'lat': ds_satellite_grid['lat'].values,
                       'lon': ds_satellite_grid['lon'].values})


# In[ ]:


# Create regridder


# In[31]:


regridder = xe.Regridder(grid_in, grid_out, 
                         method = 'bilinear', 
                         reuse_weights = True, 
                         filename = regrid_weights_dir + 'gpm_to_gfdl_cm4.nc')

# regridder = xe.Regridder(grid_in, grid_out, 
#                          method = 'bilinear', 
#                          reuse_weights = True, 
#                          filename = regrid_weights_dir + 'gpm_to_grid1.nc')


# In[32]:


dr_out = regridder(ds_satellite_grid.transpose('time','lat','lon'))


# In[13]:


# dr_out


# In[11]:


# dr_out['precipitationCal'].isel(time = 0).plot.imshow(vmin = 0, vmax = 1)


# In[10]:


# ds_satellite_grid['precipitationCal'].isel(time = 0).plot.imshow(vmin = 0, vmax = 1)


# In[12]:


# ds_model_grid['pr'].isel(time=0).plot.imshow()


# In[96]:


# dr_out['precipitationCal'].mean()


# In[ ]:





# In[ ]:


regridder = xe.Regridder(grid_in, grid_out, 
                     method = 'bilinear', 
                     reuse_weights = True, 
                     filename = regrid_weights_dir + 'gpm_to_gfdl_cm4.nc')

# regridder = xe.Regridder(grid_in, grid_out, 
#                      method = 'bilinear', 
#                      reuse_weights = True, 
#                      filename = regrid_weights_dir + 'gpm_to_grid1.nc')

for result_i in result:
    try:
        out_fname = '/'.join(result_i.split('/')[-5:])
        out_file_path = regridded_output_dir + out_fname
        rel_output_path = os.path.dirname(out_file_path)
        
        if not os.path.exists(rel_output_path):
            os.makedirs(rel_output_path)
            print(rel_output_path)
            gc.collect()
#             time.sleep(20)
        
        if not os.path.exists(out_file_path):

            ds_in = xr.open_dataset(result_i)
            ds_in['precipitationCal'] = ds_in['precipitationCal'].where(ds_in['precipitationCal'] >0, 0)
            grid_in = xr.Dataset({'lat': ds_in['lat'].values,
                               'lon': ds_in['lon'].values})


            ds_out = regridder(ds_in.transpose('time','lat','lon'))

            ds_out.to_netcdf(out_file_path)    
            
            del ds_out, ds_in

    except: 
        logging.info(f'Could not regrid {out_file_path}')
    


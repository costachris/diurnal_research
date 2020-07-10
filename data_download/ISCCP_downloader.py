#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
from ftplib import FTP, FTP_TLS
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import h5py
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

from datetime import datetime, time, timedelta
from joblib import Parallel, delayed
from mpl_toolkits.basemap import Basemap


# In[ ]:


output_data_dir = '/export/data1/cchristo/isccp_data/'
wget_url = 'https://www.ncei.noaa.gov/data/international-satellite-cloud-climate-project-isccp-h-series-data/access/isccp-basic/hgg/'


# In[ ]:


# date_range = pd.date_range(start = '2000-01',
#                           end = '2010-01', freq = '3h')
date_range = pd.date_range(start = '2010-01',
                          end = '2020-06', freq = '3h')


# In[7]:


for date_i in date_range:
    print(date_i)
    local_rel_path = output_data_dir + '%04d'%date_i.year + '/' +  '%02d'%date_i.month + '/' '%02d'%date_i.day
    if not os.path.exists(local_rel_path):
        os.makedirs(local_rel_path)
            
    wget_prefix = wget_url + '%04d'%date_i.year + '%02d'%date_i.month
    wget_filename = f'ISCCP-Basic.HGG.v01r00.GLOBAL.%04d.%02d.%02d.%02d00.GPC.10KM.CS00.EA1.00.nc'%(date_i.year,date_i.month, date_i.day, date_i.hour)
    wget_filepath =  wget_prefix + '/' +  wget_filename
    
    local_filepath = local_rel_path + '/' + wget_filename
    if not os.path.exists(local_filepath):
        result_code = os.system('wget -nc -c --retry-connrefused --waitretry=10 --quiet -P ' + local_rel_path + '/' + ' ' + wget_filepath)
#     wget_path = 


# In[ ]:


# for model_name in model_name_to_http_filename_map:
#     print(model_name)
#     all_urls_to_download = model_name_to_http_filename_map[model_name]
#     for url_to_download in all_urls_to_download:
#         local_filepath = output_data_dir + model_name + '/' + url_to_download.split('/')[-1]
#         print(local_filepath)
#         if not os.path.exists(local_filepath):
#             print("Does not exist", local_filepath)
#             result_code = os.system('wget -nc -c --retry-connrefused --waitretry=10 --quiet -P ' + output_data_dir + model_name + '/' + ' ' + url_to_download)
#     print('Done!')


# In[5]:


# date_range
# wget_prefix
# local_rel_path


# In[6]:


# wget_filename


# In[37]:


# wget_filename


# In[44]:


# local_rel_path


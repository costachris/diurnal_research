#!/usr/bin/env python
# coding: utf-8

# In[46]:


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


# In[29]:


output_data_dir = '/export/data1/cchristo/gpm_data' # no trailing backslash!


# In[47]:


ftp = FTP('arthurhou.pps.eosdis.nasa.gov')
# ftp = FTP_TLS('jsimpsonftps.pps.eosdis.nasa.gov')
# ftp = FTP_TLS('arthurhou.pps.eosdis.nasa.gov')
ftp.login('cchristo@caltech.edu','cchristo@caltech.edu')


# In[ ]:


# ftp = FTP('arthurhou.pps.eosdis.nasa.gov')
# ftp.login('cchristo@caltech.edu','cchristo@caltech.edu')


# In[48]:


ftp.nlst()


# In[ ]:





# In[55]:


ftp.cwd('/gpmdata/2001/01/01/imerg')
ftp.nlst()
# ftp.cwd('/gpmdata/2000/05/01')


# In[24]:


# sorted(ftp.nlst())
# ['3B-HHR.MS.MRG.3IMERG.20010101-S000000-E002959.0000.V06B.HDF5',


# In[ ]:


# download all files in a daterange 


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


# In[ ]:





# In[31]:


# dates = pd.date_range(start = '2000-06-01', end = '2011-06-01')
# dates = pd.date_range(start = '2000-06-01', end = '2011-06-01')
dates = pd.date_range(start = '2011-06-01', end = '2020-06-01')


# In[42]:


# dates[0].strftime('%Y/%M/%d')
retry_connection_counter = 0
for date_i in dates:
#     ftp.cwd('/gpmdata/2009/01/01/imerg')
    try:
        rel_dir_to_files = date_i.strftime('/gpmdata/%Y/%m/%d/imerg')
        logging.info(f'Downloading: {rel_dir_to_files}')
        try:
            ftp.cwd(rel_dir_to_files)
        except: 
            logging.info(f'Error with {rel_dir_to_files}!')


        sorted_fnames = sorted(ftp.nlst())

        out_file_rel_path = output_data_dir + rel_dir_to_files

        # if dir for month/day doesn't exist, create it 
        if not os.path.exists(out_file_rel_path):
            os.makedirs(out_file_rel_path)

            for fname in sorted_fnames:
                try:
                    out_file_path = out_file_rel_path + '/' + fname
                    if not os.path.exists(out_file_path):
                        with open(out_file_path, 'wb') as out_file: # Open a local file to store the downloaded file
                            print(fname)
                            ftp.retrbinary('RETR ' + fname, out_file.write, 8*1024) # Enter the filename to download
                        # open and convert to NETCDF, then save again
                        hf = h5py.File(out_file_path, 'r')
                        raw_data = hf['Grid']
                        raw_data_np = raw_data.get('precipitationCal')[()]
                        lat_np = raw_data.get('lat')[()]
                        lon_np = raw_data.get('lon')[()]
                        time_np = raw_data.get('time')[()]
                        time_np = np.array(datetime.utcfromtimestamp(time_np.item())).reshape([1])

                        # make numpy arrays into dataarray 
                        da = xr.DataArray(raw_data_np.T, dims = ('lat', 'lon', 'time'), 
                              coords = {'time': time_np,
                                        'lon': lon_np,
                                        'lat': lat_np})
                        da.name = 'precipitationCal'

                        # save as netcdf
                        new_out_file_path = out_file_rel_path + '/'+ fname[:-5] + '.nc'
                        da.to_netcdf(new_out_file_path)

                        # delete hdf5 file if it exists 
                        if os.path.exists(out_file_path):
                            os.remove(out_file_path)
                except Exception as e: 
                    logging.info(f'Failed on : {out_file_path}')
                    logging.info(f'Reason : {str(e)}')
                
            
    except Exception as f:
        logging.info(f'Failed on date : {str(date_i)}')
        logging.info(f'Reason : {str(f)}')
        
        retry_connection_counter += 1
        if retry_connection_counter < 10: 
            ftp = FTP('arthurhou.pps.eosdis.nasa.gov')
            ftp.login('cchristo@caltech.edu','cchristo@caltech.edu')
            

    


# In[37]:


# ftp.pwd()
# out_file_path.split('/')[:-1]


# In[45]:


out_file_path


# In[ ]:


# ftp = FTP('arthurhou.pps.eosdis.nasa.gov')
# ftp.login('cchristo@caltech.edu','cchristo@caltech.edu')


# In[44]:


# !wget -r "ftp://cchristo@caltech.edu:cchristo@caltech.edu@arthurhou.pps.eosdis.nasa.gov/sm/730/gpmdata/2000/06/01/imerg"


out_file_rel_path


# In[56]:


# listing = []
# ftp.retrlines("LIST", listing.append)
# words = listing[0].split(None, 8)
# filename = words[-1].lstrip()
# filename = '/sm/730/gpmdata/2000/06/01/imerg/3B-HHR.MS.MRG.3IMERG.20000601-S000000-E002959.0000.V06B.HDF5'


# In[60]:


# filename
# words
# listing
# words


# In[39]:


# listing[0]


# In[ ]:





# In[ ]:


# wget --user=cchristo@caltech.edu --password=cchristo@caltech.edu ftp://arthurhou.pps.eosdis.nasa.gov/sm/730/gpmdata/2000/06/01/imerg/3B-HHR.MS.MRG.3IMERG.20000601-S000000-E002959.0000.V06B.HDF5


# In[62]:


# from nasadap import agg

# ###############################
# ### Parameters

# cache_dir = 'nasa/cache/nz'
# save_dir = 'nasa/precip'

# username = 'cos14' # Need to change!
# password = 'Clouds14' # Need to change!

# mission = 'gpm'
# freq = 'M'
# product = '3IMERGHH'
# version = 6
# datasets = ['precipitationCal']

# min_lat=-49
# max_lat=-33
# min_lon=165
# max_lon=180
# dl_sim_count = 50
# tz_hour_gmt = 12

# agg.time_combine(mission, product, version, datasets, save_dir, username, password,
#                   cache_dir, tz_hour_gmt, freq, min_lat, max_lat, min_lon,
#                   max_lon, dl_sim_count)


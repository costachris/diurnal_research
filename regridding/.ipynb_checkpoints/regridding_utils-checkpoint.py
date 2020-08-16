import os
from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import xesmf as xe
import logging
import grid_definitions


# regrid_weights_dir = '/export/data1/cchristo/regridding_weights/'

####### General regridding utilities ########
def _regrid_generate_fname(ds_in, ds_out,
                          method, label= ''):
    '''Given descriptions of input grid and output grid, generate filename.'''
    res_in_lat = '{val:2.1f}'.format(val = np.median(np.diff(ds_in['lat'].values)))
    res_in_lon = '{val:2.1f}'.format(val = np.median(np.diff(ds_in['lon'].values)))
    res_out_lat = '{val:2.1f}'.format(val = np.median(np.diff(ds_out['lat'].values)))
    res_out_lon = '{val:2.1f}'.format(val = np.median(np.diff(ds_out['lon'].values)))
    
    in_shape_str = str(ds_in['lat'].size) + 'x' + str(ds_in['lon'].size)
    out_shape_str = str(ds_out['lat'].size) + 'x' + str(ds_out['lon'].size)
    
    
    out_string = 'res_' + res_in_lat + 'x' + res_in_lon + '_to_' + res_out_lat + 'x' + res_out_lon + \
        '_size_' + in_shape_str + '_to_' + out_shape_str + '_' + method + '_' + label + '.nc'
    return out_string
    
    
def regrid(ds_in, 
           ds_out,
           method = 'bilinear', 
           label = '',
           reuse_weights = True,
           field_name = None, 
           regrid_weights_dir = None,
           regridder = None):
    '''Regrid input dataset to output dataset. 
    If regridding weights already exist in regrid_weights_dir, use them. 
    Otherwise create new ones.
    Note: only works for rectilinear grids!  Input lat/lon coords should be 1D.
    
    Args
    ---------
        ds_in - ds.DataSet
            xarray dataset that contains coordinates ('lon', 'lat')
            Input lat/lon coords should be 1D
        reuse_weights - bool 
            If weights already exist, use them to regrid. If not, create, save, then use.
        regrid_weights_dir
            path to search for weight. Don't forget trailing backslash!
        field_name - str
            Field to regrid. If `None`, all fields will be regridded. 
    
    '''
    if (ds_in['lat'].ndim > 1) | (ds_in['lon'].ndim > 1):
        raise(Exception('Input lat/lons are not 1d!'))
    
    if (ds_out['lat'].ndim > 1) | (ds_out['lon'].ndim > 1):
        raise(Exception('Output lat/lons are not 1d!'))
    
    out_weights_filepath = regrid_weights_dir + \
                            _regrid_generate_fname(ds_in, ds_out, method = method, label = label)
    
    if not (field_name is None):
        ds_in = ds_in[field_name].to_dataset()
#     print(out_weights_filepath)
    if reuse_weights:
        if not regrid_weights_dir:
            raise(Exception('regrid_weights_dir not specified!'))
        # if weight file exists, use it
        if os.path.exists(out_weights_filepath):
            regridder = xe.Regridder(ds_in, ds_out, 
                         method = method, 
                         reuse_weights = True, 
                         filename = out_weights_filepath)
            ds_regridded = regridder(ds_in)
            return ds_regridded
        # if weight file doesn't exist, create it
        else:
            print('Weight file does not exist. Computing...')
#             rel_weight_path = os.path.dirname(out_weights_filepath)
            if not os.path.exists(regrid_weights_dir):
                os.makedirs(regrid_weights_dir)
            
            regridder = xe.Regridder(ds_in, ds_out, 
                         method = method, 
                         filename = out_weights_filepath)
            
            ds_regridded = regridder(ds_in)
            return ds_regridded
        
#     if (not reuse_weights) & (not regridder):
#         regridder = xe.Regridder(ds_in, ds_out, 
#                              method = method, 
#                              reuse_weights = False)




###### Utilities for regridding NetCDFs in a directory

def regrid_cmip_models_in_dir(model_list, 
                              ds_out, 
                              input_data_dir_rel_path, 
                              output_data_dir_rel_path,
                              regrid_weights_dir, 
                              fname = '1985-01_2006-01_precip.nc',
                              method = 'nearest_s2d',
                              out_grid_name = 'grid1',
                              field_name = None,
                              overwrite = True):
    '''Regrids netCDFs with name fname in subdirectories (model names in model_list)
    Args
    -----
    model_list - list of str
        list containing cmip model names
    out_grid_name - str
        used to label subdirs, files in output
    
    '''
    for model_ii in model_list:
        in_ds_path = input_data_dir_rel_path + model_ii + '/' + fname 
        try:
            ds_in = xr.open_dataset(in_ds_path)
            ds_regridded = regrid(ds_in, 
                                  ds_out,
                                  method = method,
                                  label = model_ii + '_to_' + 'grid1',
                                  regrid_weights_dir = regrid_weights_dir,
                                  field_name = field_name)
            out_ds_rel_path = output_data_dir_rel_path + out_grid_name + '/' + model_ii + '/'
            if not os.path.exists(out_ds_rel_path):
                os.makedirs(out_ds_rel_path)
            ds_regridded.to_netcdf(out_ds_rel_path + out_grid_name  + '_' + fname)

        except Exception as e:
            print('Could not regrid ' + model_ii + ' : ' + str(e))
'''Functions for computing stats on diurnal cycle indicies after being computed'''

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import circmean, circstd, mode
import random as rd
from mpl_toolkits.basemap import Basemap
import taylorDiagram
from phaseDiagram import *

import diurnal_config

from diurnal_utils import *
from fetch_model_helper import *

# %run cmip_metrics.py

import matplotlib 

FLUX_TO_MM_HR = 60*60
HOURS_TO_RADIANS = 2*np.pi/24

FLUX_TO_MM_HR = 60*60
MM_HR_TO_MM_DAY = 24.0
MM_HR_TO_MM_YR = 24*365
FLUX_TO_MM_YR = FLUX_TO_MM_HR*MM_HR_TO_MM_YR
FLUX_TO_MM_DAY = FLUX_TO_MM_HR*MM_HR_TO_MM_DAY


def circdiff(A,B, max_val = 24.0):
    '''Calculate absolute difference between two circular quantities '''
#     A_rad, B_rad = A*HOURS_TO_RADIANS, B*HOURS_TO_RADIANS 
    abs_diff =  max_val/2 - abs(abs(A - B) - max_val/2)
    return abs_diff

def circrmse(A, B, max_val = 24.0):
    diff = circdiff(A, B, max_val = 24.0)
    return np.sqrt((np.sum((diff)**2)/len(diff)))

def phase_circmean(arr):
    '''Use scipy circmean to calcualte circular mean of array of hours [0, 24]'''
    return circmean(arr, low =  0.0, high = 24.0)


def phase_circmean_internal(samples, high=24.0, low=0.0):
    '''Compute mean of circular quantity. '''
    samples = np.asarray(samples)
    
    sin_samp = np.sin((samples - low)*2.*np.pi / (high - low))
    cos_samp = np.cos((samples - low)*2.*np.pi / (high - low))
    
#     samples, sin_samp, cos_samp, nmask = _circfuncs_common(samples, high, low,
#                                                            nan_policy=nan_policy)

    sin_sum = sin_samp.sum()
    cos_sum = cos_samp.sum()
    res = np.arctan2(sin_sum, cos_sum)
    
    # if angle is <0 (clockwise on unit circle), 
    # switch to equivalent positive (counter clockwise) angle
    if res < 0:
        res = 2*np.pi + res
    
    return res*(high - low)/2.0/np.pi + low
    


def lat_weighted_circ_mean(df, 
                           field_name, 
#                            weight_field_2 = 'pr_mean', 
                           weight_field_2 = None, 
                           high=24.0, 
                           low=0.0):
    '''Given dataframe with circular quantity, find lat weighted mean.
    weight_field_2.
    
    Args
    -----------
    
    weight_field_2 - str
        Field name of additional field to weight mean by.
    
    Returns
        Weighted mean of field 
    '''
    samples = np.asarray(df[field_name].values)
    
    sin_samp = np.sin((samples - low)*2.*np.pi / (high - low))
    cos_samp = np.cos((samples - low)*2.*np.pi / (high - low))
    
    if 'lat' in df.index.names:
        df = df.reset_index()
    
    weights = np.cos(np.deg2rad(df['lat']))
    
    if weight_field_2:
        print('double-weighted')
        weights = np.multiply(df[weight_field_2].values, weights)
        
#     sin_ave = np.average(sin_samp, weights = weights)
#     cos_ave = np.average(cos_samp, weights = weights)
#     res = np.arctan2(sin_ave, cos_ave)

    weights_sum = weights.sum()
    
    num_sum_sin = (weights * sin_samp).sum()
    num_sum_cos = (weights * cos_samp).sum()
    
    res = np.arctan2(num_sum_sin/weights_sum, num_sum_cos/weights_sum)
    
    # if angle is <0 (clockwise on unit circle), 
    # switch to equivalent positive (counter clockwise) angle
    if res < 0:
        res = 2*np.pi + res
    
    return res*(high - low)/2.0/np.pi + low
    

def hour_circstd(arr):
    '''Use scipy circstd to calcualte circular mean of array of hours [0, 24]'''
    return circstd(arr, low = 0.0, high = 24.0)

def sin_hour(hours):
    '''Give array of hours, convert sin quanitity to deal with cirular quantity'''
    return np.sin(hours*HOURS_TO_RADIANS)

def hour_corr(arr1, arr2):
    '''Find circular correlation between two arrays of hours [0, 24]'''
    return np.corrcoef(sin_hour(arr1), sin_hour(arr2))[0,1]

def ampl_weighted_mean_func(df):
    return (df['precip_weights']*df['ampl_season']).sum()

def mode_apply(df):
    '''Scipy mode returns both mode and counts, os wrap function to only return mode. '''
    return mode(df)[0].item()

def df_mean_lat_weighted(df, field_name):
    '''Calculate mean of field in Dataframe, weighting by cos(latitude)
    
    
    Args
    -----
    df - pd.DataFrame
        Dataframe containing field to find mean of (with id = `field_name`) and `lat` 
        as an index or standard column. 
        
    Returns
    --------
    weighted_mean - float
    
    '''
    if 'lat' in df.index.names:
        df = df.reset_index()
    
    weights = np.cos(np.deg2rad(df['lat']))
    
    num_sum = (weights * df[field_name]).sum()
    denom_sum = weights.sum()
    
    return num_sum/denom_sum


# def mode_apply_phase(df, round_precision = 1):
#     '''Handle binning of data before taking mode for given variables'''
#         df.round(round_precision)
# def mode_apply_ampl(df, round_precision = 4):

def _open_and_preprocess_gpm(input_data_dir_gpm,
                             ds_land_sea,
                             yearly_mean_bool):
    ds_gpm = xr.open_dataset(input_data_dir_gpm +  'grid1_2000-06_2016-06_precip.nc')

    if yearly_mean_bool:
        ### compute means
        # # take mean across seasons using circular mean for phase
        ds_gpm_phase_year_mean = xr.apply_ufunc(phase_circmean, ds_gpm['phase_season'], 
        #                                    kwargs = {'low' : 0.0, 'high' : 24.0},
                                           input_core_dims=[["season"]], 
                                           vectorize = True)
        # compute yearly mean other data 
        ds_gpm = ds_gpm.mean(dim = 'season')
        ds_gpm['phase_season'] = ds_gpm_phase_year_mean

    #########

    ds_gpm['land_sea_mask'] = ds_land_sea['GLDAS_mask']
    df_gpm = ds_gpm.to_dataframe()
    
    return df_gpm


def _load_all_means(cmip_rel_dir,
                    filename = 'grid1_1985-01_2006-01_mean.nc'):
    '''Given path to mean precip fields, where models are subdirs, load and return a dataframe'''
    model_names = os.listdir(cmip_rel_dir)
    mean_dfs = []
    for model_name in model_names: 
        ds_i = xr.open_dataset(cmip_rel_dir + model_name + '/' + filename)
        ds_i_df = ds_i.to_dataframe()
        if 'bnds' in ds_i_df.index.names:
            ds_i_df = ds_i_df.reset_index('bnds').drop(['lat_bnds', 'lon_bnds', 'bnds'], axis = 1)


        ds_i_df['model_name'] = model_name
        ds_i_df = ds_i_df.rename({'pr':'pr_mean'}, axis = 1)
#         ds_i_df = ds_i_df.reset_index()
#         ds_i_df = ds_i_df.reset_index().set_index(['model_name','lat','lon'])
        mean_dfs.append(ds_i_df)
    
    return pd.concat(mean_dfs, axis = 0)

    
    

### filtering functions for pd.DataFrames 
def filter_by_lat(df, min_lat, max_lat, absolute_value = False):
    '''Filter pd.DataFrame by lat. 
    Args
    -------
    
    absolute_value - bool
        use absolute value of min/max lat to filter. Ie. if
        min_lat = 35 and max_lat = 60, find regions in both northern
        and southern hemisphere. (35 to 60 & -35 to -60). 
    '''
    if 'season' in df.index.names:
        desired_index_names = ['season','lat','lon']
        df = df.reset_index(desired_index_names)
        
    else:
        desired_index_names = ['lat','lon']
        df = df.reset_index(desired_index_names)
    
    if absolute_value:
        return df[(abs(df['lat'])>= min_lat) & (abs(df['lat']) <= max_lat)].set_index(desired_index_names)

    else:
        return df[(df['lat']>= min_lat) & (df['lat'] <= max_lat)].set_index(desired_index_names)

def filter_by_season(df, season):
    '''Filter df by season. '''
    df = df.reset_index('season')
    return df[df['season'] == season]
    
    
####### ds operators 
def _create_mask_for_all_models(model_list,
                                **kwargs):
    var_dict = {}
    for model_name in model_list: 
        var_dict[model_name] = _create_mask_from_model(model_name = model_name,
                                            **kwargs)
    out_df = pd.concat(var_dict)
    # label model index and reset
    idx_names = list(out_df.index.names)
    idx_names[0] = 'model_name'
    out_df.index.names = idx_names
    
    out_df = out_df.reset_index(['model_name'])
    
    return out_df

def _create_mask_from_model(data_dir, 
                            model_name,
                            field_to_mask = 'pr',
                            field_threshold = 1.0,
                            unit_conversion = FLUX_TO_MM_DAY,
                            operator = 'lt', 
                            file_name = 'grid1_1985-01_2006-01_mean.nc'):
    '''
    Given path to model NetCDF, field, and theshold. Create mask from model.
    Args
    ------
    data_dir - str
        dir containing cmip models in seperate folders
    model_name - str
        Name of model to compute mask for
    field_to_mask - str
        Name of field to compute mask
    field_threshold - float
        Field threshold
    operator - str {'lt', 'gt'}
        Mask based on values < ('lt') or > ('gt') threshold
        
    Returns
    --------
        mask - pd.DataFrame 
            Mask containing 1 where condition is true and 0 elsewhere. 
    '''
    
    ds_field = xr.open_dataset(data_dir + model_name + '/' + file_name)
    
    if not (unit_conversion is None):
        ds_field[field_to_mask] = ds_field[field_to_mask]*unit_conversion
    
    # mask 
    if operator == 'lt':
        field_where = ds_field[field_to_mask].where(ds_field[field_to_mask] < field_threshold)
    elif operator == 'gt':
        field_where = ds_field[field_to_mask].where(ds_field[field_to_mask] > field_threshold)


    field_mask = (field_where > 0)
    field_mask.name = 'mask'
    field_mask_df = field_mask.to_dataframe()
   
    return field_mask_df

def _get_mean_field(data_dir,
                    field_name = 'rlut',
                    file_name = '1985-01_2006-01_mean.nc',
                    filter_lat = False,
                    landsea_mask_df = None,
                    var_mask_df = None,
                    landsea_bool = None,
                    weighted_mean_bool = True,
                    **kwargs):
    '''Given path to cmip files (with subdirs for each model),
    take mean of each model and return dataframe.
    Args
    -----
    data_dir - str
        dir containing cmip models in seperate folders
        
    filter_lat - bool
        whether or not to filter by lat. 
    landsea_mask_df - pd.DataFrame
        Pandas dataframe containing lat/lon as multiindex keys and 
        mask as column (0 for water, 1 for land). df should only
        contain one column. 
    var_mask_df - pd.DataFrame or 'auto'
        only compute averages where mask is 1. Either provide mask 
        dataframe or 'auto'. Setting 'auto' will use
        `_create_mask_from_model` to automatically generate a seperate
        mask for each model. 
    landsea_bool - bool 
        If 1 return land only, if 0 return sea only. 
    weighted_mean_bool - bool
        Weight field by cos(lat) if true when computing mean. 
            
    **kwargs
        keyword arguments to be past to lat filtering function
        (filter_by_lat)
    
    '''
    var_dict = {}

    model_names = os.listdir(data_dir)
#     print(model_names)
    for model_name in model_names:
        try:
            ds_field = xr.open_dataset(data_dir + model_name + '/' + file_name)
            df_field = ds_field.to_dataframe()
            
            ## if landsea mask filter
            if (not (landsea_mask_df is None)) & \
               (not (landsea_bool is None)): 
                df_field_mask = pd.merge(df_field, landsea_mask_df,
                                         how = 'left',
                                         left_index = True, 
                                         right_index = True)
                landsea_mask_id = landsea_mask_df.columns[0]
                
                # select land
                if landsea_bool == 1:
                    df_field_mask = \
                        df_field_mask[df_field_mask[landsea_mask_id] == 1]
                
                # select water
                if landsea_bool == 0:
                    df_field_mask = \
                        df_field_mask[df_field_mask[landsea_mask_id] == 0]
                
                
                 # apply other mask, if specified
                if not (var_mask_df is None):
                    if type(var_mask_df) == dict:
                        var_mask_df = _create_mask_from_model(
                                        model_name = model_name,
                                        **var_mask_df)
                        
                    var_mask_id = var_mask_df.columns[0]
                    df_field_mask = pd.merge(df_field_mask, var_mask_df,
                                         how = 'left',
                                         left_index = True, 
                                         right_index = True)
                    
                    df_field_mask = \
                        df_field_mask[df_field_mask[var_mask_id] == 1]
                
                    
                df_field = df_field_mask.loc[:,field_name]
                
               
                
                
            if filter_lat:
                df_field_filtered = filter_by_lat(df_field, **kwargs)
                if weighted_mean_bool:
                    var_dict[model_name] = df_mean_lat_weighted(df_field_filtered, field_name)
                else:
                    var_dict[model_name] = df_field_filtered[field_name].mean()
#                 return df_field_filtered
            
            else:
                if weighed_mean_bool: 
                    var_dict[model_name] = df_mean_lat_weighted(df_field, field_name)
                else:
                    var_dict[model_name] = df_field[field_name].mean()
            

        except Exception as e: 
            print('Could not process ', model_name, ' ', str(e))
    out_df = pd.DataFrame.from_dict(var_dict, orient = 'index')
    out_df.columns = [field_name,] 
    return out_df

def _get_mean_fields(
    field_names,
    mean_fields_to_rel_path_map,
    var_mask_df = None,
    **kwargs):
    '''Apply _get_mean_field over a sequence of variables
    and return a single DataFrame.
    
    Args
    -------
    rel_dir - str
        Relative directory to files containing mean fiels 
    '''
    out_df = pd.DataFrame()
    dfs_to_concat = []
    for field_name in field_names: 
        data_dir = mean_fields_to_rel_path_map[field_name]
        # open dataset with fields to 
        dfs_to_concat.append(_get_mean_field(
                            data_dir = data_dir, 
                            field_name = field_name, 
                            var_mask_df = var_mask_df,
                            **kwargs))
    return pd.concat(dfs_to_concat, axis = 1)


    
def _merge_models_into_df(model_names, input_data_dir,
                          filename = '1985-01_2006-01_precip.nc',
                          compute_year_mean_bool = False, 
                          verbose = False):
    '''Given list of model names and path to cmip files containing seasonal components, 
    open and take yearly mean (ensuring a circular mean is used for phase).'''
    objs = {}
    df = pd.DataFrame()
    for model_name_i in range(len(model_names)):  
        try:
            model_name = model_names[model_name_i]
            path_to_file = input_data_dir + model_name + '/' + filename
            ds_i = xr.open_dataset(path_to_file)
            if compute_year_mean_bool:
                ds_i_phase_year_mean = xr.apply_ufunc(phase_circmean, ds_i['phase_season'], 
#                                    kwargs = {'low' : 0.0, 'high' : 24.0},
                                   input_core_dims=[["season"]], 
                                   vectorize = True,
                                   dask = 'allowed')
                ds_i = ds_i.mean(dim = 'season')
                ds_i['phase_season'] = ds_i_phase_year_mean
            objs[model_names[model_name_i]] = ds_i.to_dataframe()
            objs[model_names[model_name_i]]['model_name'] = model_name
        
        except Exception as e:
            if verbose:
                print('Could not get ', model_name_i, ' : ', e)
    df = pd.concat(list(objs.values()), axis = 0)
    return df


def _merge_df_on_model_name(df1, df2):
    '''Merge Dataframes on lat, lon, and model name. Assumes 
    ('lat', 'lon') are Multilevel indices in both Dataframes.'''
    df1 = df1.reset_index()
    
    df_merged = pd.merge(df1, df2,
                         how = 'left',
                         on = ['lat', 'lon', 'model_name'])
    return df_merged.set_index(['lat','lon'])
    
    
    
def full_analysis(df_for_stats,
                  df_for_stats_true,
                  cmip_sensitivities,
                  field_means_df = None,
                  field = 'phase_season',
                  agg_method = 'mode',
                  error_stats = False,
                  var_mask_df = None,
                  min_lat = None, 
                  max_lat = None, 
                  absolute_value = False,
                  ):
    
#     if not (season is None):
#         df_for_stats = filter_by_season(df_for_stats, season)
#         df_for_stats_true =  filter_by_season(df_for_stats_true, season)
    
    df_for_stats = df_for_stats.copy()
    df_for_stats_true = df_for_stats_true.copy()
    
    df_for_stats.loc[:,'ampl_season'] = df_for_stats['ampl_season'].apply(lambda x: x*FLUX_TO_MM_DAY)
    df_for_stats.loc[:,'mu_season'] = df_for_stats['mu_season'].apply(lambda x: x*FLUX_TO_MM_DAY)
    
    # split obs into land/sea 
    df_for_stats_true_land = df_for_stats_true[df_for_stats_true['land_sea_mask'] == 1]
    df_for_stats_true_water = df_for_stats_true[df_for_stats_true['land_sea_mask'] == 0]

    # split model into land/sea
    df_for_stats_land = df_for_stats[df_for_stats['land_sea_mask'] == 1]
    df_for_stats_water = df_for_stats[df_for_stats['land_sea_mask'] == 0]
    
    
    if (not min_lat is None) & (not max_lat is None):
        df_for_stats_land = filter_by_lat(df_for_stats_land, 
                                          min_lat, 
                                          max_lat, 
                                          absolute_value=absolute_value)
        
        df_for_stats_true_land = filter_by_lat(df_for_stats_true_land, 
                                               min_lat, 
                                               max_lat, 
                                               absolute_value=absolute_value)


        df_for_stats_water = filter_by_lat(df_for_stats_water, 
                                           min_lat, 
                                           max_lat, 
                                           absolute_value=absolute_value)
        
        df_for_stats_true_water = filter_by_lat(df_for_stats_true_water, 
                                                min_lat, 
                                                max_lat, 
                                                absolute_value=absolute_value)
#     return df_for_stats_water, None
    # apply other mask, if specified
    if not (var_mask_df is None):
        
#         var_mask_id = var_mask_df.columns[0]
        var_mask_id = 'mask'
        # merge seperate mask for each model
        if 'model_name' in var_mask_df:
            df_for_stats_water = _merge_df_on_model_name(df_for_stats_water, var_mask_df)
        # TODO: if different mask for each model, which to use for satellite?
#             df_for_stats_true_water = _merge_df_on_model_name(df_for_stats_true_water, var_mask_df)
        
        # merge single mask for all models
        else:
            df_for_stats_water = pd.merge(df_for_stats_water, var_mask_df,
                                 how = 'left',
                                 left_index = True, 
                                 right_index = True)

#             df_for_stats_true_water = pd.merge(df_for_stats_true_water, var_mask_df,
#                                  how = 'left',
#                                  left_index = True, 
#                                  right_index = True)

        df_for_stats_water = \
            df_for_stats_water[df_for_stats_water[var_mask_id] == 1]
        
#         df_for_stats_true_water = \
#             df_for_stats_true_water[df_for_stats_true_water[var_mask_id] == 1]

#     return df_for_stats_water, None
    

    model_error_stats_df_land = compute_stats(df_for_stats_land,
                 df_for_stats_true_land,
                 field = field,
                 agg_method = agg_method,
                 error_stats = error_stats,                           
                 additional_stats = True,)

    ### compute stats for land/water
    model_error_stats_df_water = compute_stats(df_for_stats_water,
                     df_for_stats_true_water,
                     agg_method = agg_method,
                     error_stats = error_stats,
                     field = field,
                     additional_stats = True,)
#     return(model_error_stats_df_water, cmip_sensitivities)
    if 'season' in model_error_stats_df_water.index.names:
        model_error_stats_df_water = model_error_stats_df_water.reset_index('season')
        
    if 'season' in model_error_stats_df_land.index.names:
        model_error_stats_df_land = model_error_stats_df_land.reset_index('season')
        
    all_stats_df_water = pd.merge(model_error_stats_df_water, cmip_sensitivities, 
                             how = 'left',
                             left_index = True, 
                             right_index = True)
    all_stats_df_land = pd.merge(model_error_stats_df_land, cmip_sensitivities, 
                         how = 'left',
                         left_index = True, 
                         right_index = True)

    if not field_means_df is None:
        all_stats_df_water = pd.merge(all_stats_df_water, field_means_df,
             how = 'left',
             left_index = True, 
             right_index = True)
        all_stats_df_land = pd.merge(all_stats_df_land, field_means_df,
             how = 'left',
             left_index = True, 
             right_index = True)
            
    return (all_stats_df_water, all_stats_df_land)

    
# compute stats on df containing indicies 

def compute_stats(df_for_stats,
                 df_for_stats_true,
                 field = 'phase_season',
                 agg_method = 'mean',
                 error_stats = True,
                 additional_stats = True, 
                 ecs_dict = None,
                 tcr_dict = None,
                 rlut_dict = None,
                 rsut_dict = None,
                 pr_dict = None,
                 clt_dict = None):
    '''Given pd.DataFrames containing cmip model and satellite diurnal cycle indicies, compute various error statistics
    and return dataframe. Circular statistics are used for the phase field. Note: when agg method = `mean`
    quantities are weighted by cos(lat)
    
    Args
    -----------
        df_for_stats - pd.DataFrame
            Dataframe containing CMIP model diurnal stats
        df_for_stats_true - pd.DataFrame
            Dataframe containing data to validate agaisnt. 
        
        field - str
        
        error_stats - bool
            Whether or not to compute accuracy metrics (rmse, corr)
            agaisnt satellite observations. 
        
        agg_method - str in {'mean', 'mode'}
            method for aggregating field. Means are circular for phase. 
            When agg method = `mean` quantities are weighted by cos(lat).
        
        additional_stats - bool
            #TODO: right now this has to be true!
            whether to compute additional stats for comparing against ECS, TCR
        ecs_dict - dict
            Dictionary mapping model name to ECS. 
        
        
    Returns
    ----------
        model_error_stats_df - pd.DataFrame
            Dataframe containing all error stats
    
    '''
    df_for_stats = df_for_stats.copy()
    df_for_stats_true = df_for_stats_true.copy()
    
    model_error_stats = {}
    
    # logic for handling seasonal breakdowns
    if ('season' in df_for_stats.index.names) | \
       ('season' in df_for_stats.columns):
        groupby_vars = ['season', 'model_name']
        season_bool = True
    else:
        groupby_vars = 'model_name'
        season_bool = False
    
    
    if 'season' in df_for_stats.index.names:
        df_for_stats = df_for_stats.reset_index('season')
    if 'season' in df_for_stats_true.index.names:
        df_for_stats_true = df_for_stats_true.reset_index('season')
        
    # subset columns before doing pandas aggs to speed up calculations
    if 'pr_mean' in df_for_stats.columns:
        _variables_of_interest = ['ampl_season', 'phase_season', 'model_name', 'pr_mean']
    else: 
        _variables_of_interest = ['ampl_season', 'phase_season', 'model_name']
        
    if 'season' in df_for_stats.columns:
        _variables_of_interest.append('season')
    
    phase_bin_precision = 1
    ampl_bin_precision = 4
    
    if agg_method == 'mean':

#         agg_field_by_model = df_for_stats.groupby(groupby_vars).mean()
        agg_field_by_model = df_for_stats.groupby(groupby_vars).apply(df_mean_lat_weighted, field_name = 'ampl_season')
        agg_field_by_model = pd.DataFrame(agg_field_by_model)
        agg_field_by_model = agg_field_by_model.rename({0:'ampl_season'}, axis = 1)
        
        agg_phase_by_model = df_for_stats[_variables_of_interest].groupby(groupby_vars).apply(lat_weighted_circ_mean, field_name = 'phase_season')
        agg_phase_by_model = pd.DataFrame(agg_phase_by_model)
        agg_phase_by_model = agg_phase_by_model.rename({0:'phase_season'}, axis = 1)
        
        
    elif agg_method == 'mode':
#         print('here')
        agg_field_by_model = df_for_stats[_variables_of_interest].round(ampl_bin_precision).groupby(groupby_vars).agg(mode_apply)
        
        agg_phase_by_model = df_for_stats[_variables_of_interest].round(phase_bin_precision).groupby(groupby_vars).agg(mode_apply)
#     return agg_field_by_model
    if 'season' in agg_phase_by_model.index.names:
        agg_field_by_model =  agg_field_by_model.reset_index('season')
        agg_phase_by_model = agg_phase_by_model.reset_index('season')

#     return agg_field_by_model
        df_for_stats_true_grouped = df_for_stats_true.groupby('season')

    for indexers, df_i in df_for_stats.groupby(groupby_vars):

        if season_bool:
            season_i = indexers[0]
            model_name = indexers[1]
            df_for_stats_true_i = df_for_stats_true_grouped.get_group(season_i)
            df_true_field = df_for_stats_true_i[field]
        else:
            model_name = indexers
            df_true_field = df_for_stats_true[field]
        

        if error_stats:
            rmse_i = circrmse(df_i[field], df_true_field)
            model_i_corr = np.corrcoef(sin_hour(df_i[field].values), 
                                   sin_hour(df_true_field.values))[0,1]
            model_i_std = circstd(df_i[field].values, low = 0.0, high = 24.0)
        else:
            rmse_i = np.nan
            model_i_corr = np.nan
            model_i_std = np.nan

        

        

       


        # To Do
        if additional_stats == True:
         # use regular mean
#             return(agg_field_by_model, season_i)
            if 'season' in agg_field_by_model.columns:
                ampl_mean = agg_field_by_model[agg_field_by_model['season'] == season_i].loc[model_name]['ampl_season']
            else:
                ampl_mean = agg_field_by_model['ampl_season'][model_name]

            if 'season' in agg_phase_by_model.columns:
                phase_mean = agg_phase_by_model[agg_phase_by_model['season'] == season_i].loc[model_name]['phase_season']
            else:
                phase_mean = agg_phase_by_model['phase_season'][model_name]
            
        
#         model_error_stats[indexers] = [model_i_std, model_i_corr, rmse_i, 
#                                          ampl_mean, phase_mean]
        
#     model_error_stats_df = pd.DataFrame(model_error_stats).T
#     model_error_stats_df.columns = ['std', 'corr', 'rmse', 
#                                 'ampl_mean', 'phase_mean']
#     return model_error_stats
        
        
        
        
#             if (not ecs_dict is None) & (model_name in ecs_dict):
#                 ecs_i = ecs_dict[model_name]
#             else:
#                 ecs_i = np.nan
                
#             if (not tcr_dict is None) & (model_name in tcr_dict):
#                 tcr_i = tcr_dict[model_name]
#             else:
#                 tcr_i = np.nan
                
#             if (not rlut_dict is None) & (model_name in rlut_dict):
#                 rlut_i = rlut_dict[model_name]
#             else:
#                 rlut_i = np.nan
            
#             if (not rsut_dict is None) & (model_name in rsut_dict):
#                 rsut_i = rsut_dict[model_name]
#             else:
#                 rsut_i = np.nan
                
#             if (not pr_dict is None) & (model_name in pr_dict):
#                 pr_i = pr_dict[model_name]
#             else:
#                 pr_i = np.nan
            
#             if (not clt_dict is None) & (model_name in clt_dict):
#                 clt_i = clt_dict[model_name]
#             else:
#                 clt_i = np.nan

        model_error_stats[indexers] = [model_i_std, model_i_corr, rmse_i, 
                                         ampl_mean, phase_mean,] 
#     ecs_i, tcr_i,
#                                          rlut_i, rsut_i, pr_i, clt_i]

    model_error_stats_df = pd.DataFrame(model_error_stats).T
    model_error_stats_df.columns = ['std', 'corr', 'rmse', 
                                    'ampl_mean', 'phase_mean',] # 'ecs', 'tcr',
#                                     'rlut', 'rsut', 'pr', 'clt']
    
    if len(groupby_vars) == 2: 
        model_error_stats_df.index.names = ['season', 'model_name']
    return model_error_stats_df

################ Plotting #############################
def land_sea_histogram(df, 
                       nbins = 120,
                       title = None,
                       ylabel = 'Probability Density',
                       xlabel = None,
                       cmip_identifier = 'CMIP6', 
                       new_fig = True, 
                       ax = None,
                       field_id = 'phase_season'):
    '''Given a dataframe, plot histogram and corresponding pdf with land/ocean breakdown'''
    if new_fig & (not ax is None):
        plt.figure()
    
    field_id_to_xlabel = {'phase_season': 'Precipitation Phase [hours]',
                          'ampl_season': 'Amplitude'}
    if 'cmip_identifier' in df.columns:
        df_water = df[(df['land_sea_mask'] == 0) & (df['cmip_identifier'] == cmip_identifier)][field_id]
        
        df_land = df[(df['land_sea_mask'] == 1) & (df['cmip_identifier'] == cmip_identifier)][field_id]
        sns.distplot(df_water.values , label = cmip_identifier + ' ' + 'Water', bins = nbins, ax = ax)
        sns.distplot(df_land.values , label = cmip_identifier + ' ' + 'Land', bins = nbins, ax = ax)
    else: 
        sns.distplot(df[(df['land_sea_mask'] == 0)][field_id].values, 
                     label = 'Water', 
                     bins = nbins, 
                     ax = ax)
        sns.distplot(df[(df['land_sea_mask'] == 1)][field_id].values, 
                     label = 'Land', 
                     bins = nbins, 
                     ax = ax)
    if not (ax is None):
        if field_id == 'phase_season':
            ax.set_xlim([0, 24])
            ax.set_ylim([0, 0.22])
        if xlabel:
            ax.set_xlabel(field_id_to_xlabel[field_id], fontweight = 'bold')
            
        ax.set_ylabel(ylabel, fontweight = 'bold')
        if title:
            ax.set_title(title, fontweight = 'bold')
        ax.grid()
        ax.legend()
        
    else: 
        if field_id == 'phase_season':
            plt.xlim([0, 24])
            plt.ylim([0, 0.25])
        plt.xlabel(field_id_to_xlabel[field_id], fontweight = 'bold')
        plt.ylabel(ylabel, fontweight = 'bold')
        plt.title(title, fontweight = 'bold')
        plt.grid()
        plt.legend()
    
#     return plt.gca()



def plot_corr_matrix(corr_mat_ds,
                    title = 'Correlation Matrix'):
    plt.figure(figsize = (12,10))
    if title:
        plt.title(title)
    upper_tr_mask = np.triu(corr_mat_ds.corr())
    sns.heatmap(corr_mat_ds.corr(), annot = True, cmap = 'RdBu_r',
                vmin = -1, vmax = 1, center = 0, fmt='.2g',
                mask = upper_tr_mask)
    plt.tight_layout()
    

    
def summary_stats_for_df(df, 
                         agg_method = 'mode',
                         phase_bin_precision = 1,
                         ampl_bin_precision = 4):
    '''Given a pd.DataFrame, compute summary stats for 
    amplitude and phase (`mean`, `mode`, `precip_mean`).
    
    Args
    ------
    agg_method - str
        {`mean`, `mode`, `precip_mean`}
        mode : mode
        mean : cos(lat) weighted mean
        precip_mean :  cos(lat) + precip weighted mean
    
    
    '''
    if agg_method == 'mode': 
        ampl_agg = mode_apply(df['ampl_season'].round(ampl_bin_precision))
        phase_agg = mode_apply(df['phase_season'].round(phase_bin_precision))
    if agg_method == 'mean':
        ampl_agg = df_mean_lat_weighted(df, 'ampl_season')
        phase_agg = lat_weighted_circ_mean(df, 'phase_season')
        
    if agg_method == 'precip_mean':
        ampl_agg = df_mean_lat_weighted(df, 'ampl_season')
        phase_agg = lat_weighted_circ_mean(df, 'phase_season')
        
    return (ampl_agg, phase_agg)



def make_phase_plot(water_df,
                    land_df,
                    obs_water_df,
                    obs_land_df, 
                    agg_method = 'mean',
                    title = r'Diurnal Phase [hr] & Amplitude [$\frac{mm}{day}$] for CMIP6 and IMERG',
                    y_lim = (0, 1),
                    figsize = (13,8),
                    markersize = 2,
                    textsize = 8,
                    normalize_ampl = False,
                    ampl_unit_conversion_factor = 1.0):
    
    fig = plt.figure(figsize = figsize)

    
    taylor_diag = PhaseDiagram(fig = fig, 
                              label = 'IMERG', 
                              y_lim=y_lim,
                              radial_label_pos = 0
                              )
    taylor_diag.add_grid()
    
    # calculate obs phase mode
    
    ampl_observed_water, phase_observed_water = summary_stats_for_df(obs_water_df, 
                                                                     agg_method = agg_method)
    
    ampl_observed_land, phase_observed_land = summary_stats_for_df(obs_land_df, 
                                                                   agg_method = agg_method)
    
    
#     ampl_observed_water = mode_apply(obs_water_df['ampl_season'].round(4))
#     ampl_observed_land = mode_apply(obs_land_df['ampl_season'].round(4))
    
#     phase_observed_water = mode_apply(obs_water_df['phase_season'].round(1)
#     phase_observed_land = mode_apply(obs_water_df['phase_season'].round(1))

    if normalize_ampl:
        ampl_to_plot_water = 1.0
        ampl_to_plot_land = 1.0
    else:
        ampl_to_plot_water = ampl_observed_water
        ampl_to_plot_land = ampl_observed_land
    
    print(ampl_to_plot_water)
    taylor_diag.add_sample(phase = phase_observed_water, 
                               ampl = ampl_to_plot_water, 
                               marker = '*', 
                               c = 'b',
                               label = 'IMERG-Water', 
                               markersize = 13)

    taylor_diag.add_sample(phase = phase_observed_land, 
                               ampl = ampl_to_plot_land, 
                               marker = '*', 
                               c = 'g',
                               label = 'IMERG-Land', 
                               markersize = 13)


    ## Plot model phase/ampl over water
    model_list = list(water_df.index)
    for model_name_i in range(len(model_list)):

        phase_i = water_df.loc[model_list[model_name_i],:]['phase_mean']
        ampl_i = water_df.loc[model_list[model_name_i],:]['ampl_mean']

        if normalize_ampl:
            ampl_i = ampl_i/ampl_observed_water
        taylor_diag.add_sample(phase = phase_i, 
                               ampl = ampl_i, 
                               marker = None,
                               linestyle = None,
                               c = 'b',
                               label = str(model_name_i) + ': ' + model_list[model_name_i], 
                               markersize = markersize)
        taylor_diag.add_text(phase = phase_i, 
                             ampl = ampl_i,
                             text = model_name_i,
                             label = str(model_name_i) + ': ' + model_list[model_name_i], 
                             c = 'b',
                             size = textsize,
                             weight = 'bold')

    ## Plot model phase/ampl over land
    for model_name_i in range(len(model_list)):
        phase_i = land_df.loc[model_list[model_name_i],:]['phase_mean']
        ampl_i = land_df.loc[model_list[model_name_i],:]['ampl_mean']

        if normalize_ampl:
            ampl_i = ampl_i/ampl_observed_land

        taylor_diag.add_sample(phase = phase_i, 
                               ampl = ampl_i, 
                               marker = None, 
                               c = 'g',
                               linestyle = None,
    #                            label = str(model_name_i) + ': ' + model_list[model_name_i], 
                               label = None,
                               markersize = markersize)
        taylor_diag.add_text(phase = phase_i, 
                             ampl = ampl_i,
                             text = model_name_i,
                             c = 'g',
                             size = textsize,
                             weight = 'bold')

    ax = plt.gca()
    ax.text(-.1, ax.get_rmax()/2.,r'Amplitude [$\frac{mm}{day}$]',
        rotation=90, ha='center',va='center')
    plt.xlabel('Local Solar Time [Hours]', weight = 'bold')

    leg = plt.legend(loc = 'center left', bbox_to_anchor=(1.1,0.5), prop={'size': 11}, handlelength = 0, markerscale = 0.8)


    plt.title(title, weight = 'bold')
    

    
def linear_regression(X, Y):
    beta=np.linalg.solve(X.T.dot(X),X.T.dot(Y))
    return beta
def confidence_bounds(x,betas):
    predictions=[x.dot(b) for b in betas]
    predictions.sort()
    CI_95=[predictions[24], predictions[-25]]
    CI_68=[predictions[160], predictions[-160]]
    CI_99=[predictions[1], predictions[-1]]
    return CI_99,CI_95,CI_68
def regression_bounds(x,y,XN):
    X=x
    Y=y
    sorted_x = np.sort(X)
    argind = np.argsort(X)
    sorted_y = np.zeros((len(x),1))
    for i in range(len(x)):
        sorted_y[i] = Y[argind[i]]
    X=np.expand_dims(sorted_x,1)
    Y=sorted_y
    X = np.hstack((X, np.ones((X.shape[0], 1), dtype=X.dtype)))
    beta=linear_regression(X, Y)
    betas=[]
    N= X.shape[0]
    for resample in range(1000):
        index = [rd.choice(range(N)) for i in range(N)]
        re_X=X[index,:]
        betas.append(linear_regression(re_X, Y[index,:]))
    upper99=np.zeros((len(x),1))
    lower99=np.zeros((len(x),1))
    upper95=np.zeros((len(x),1))
    lower95=np.zeros((len(x),1))
    upper68=np.zeros((len(x),1))
    lower68=np.zeros((len(x),1))
    cnt=0
    for xn in XN:
        CB_99,CB_95,CB_68 = confidence_bounds(xn,betas)
        upper99[cnt]=(CB_99[1])
        lower99[cnt]=(CB_99[0])
        upper95[cnt]=(CB_95[1])
        lower95[cnt]=(CB_95[0])
        upper68[cnt]=(CB_68[1])
        lower68[cnt]=(CB_68[0])
        cnt=cnt+1
    return upper99,lower99,upper95,lower95,upper68,lower68,beta


#### Old code 


# from astropy import units as u
# model_error_stats = {}
# additional_stats = True

# mean_field_by_model = df_for_stats.groupby('model_name').mean()
# circmean_phase_by_model = df_for_stats[['phase_season', 'model_name']].groupby('model_name').apply(phase_circmean)


# HOURS_TO_RADIANS = 2*np.pi/24
# # HOURS_TO_RADIANS = 360.0/24
# for model_name, df_i in df_for_stats.groupby('model_name'):
#     diff  = df_i[field] - df_gpm[field]
#     rmse_i = circrmse(df_i[field], df_gpm[field])
    
#     model_i_corr = np.corrcoef(sin_hour(df_i[field].values), sin_hour(df_gpm[field].values))[0,1]
# #     model_i_corr = np.corrcoef(df_i[field].values, df_gpm[field].values)[0,1]
# #     model_i_corr = astro_circstats.circcorrcoef(df_i[field].values * HOURS_TO_RADIANS,  
# #                              df_gpm[field].values * HOURS_TO_RADIANS)
    
#     model_i_std = circstd(df_i[field].values, low = 0.0, high = 24.0)
    
    
#     # To Do
#     if additional_stats == True:
#         if model_name in cmip6_ecs:
#             ecs_i = cmip6_ecs[model_name]
#             # use regular mean
#             ampl_mean = mean_field_by_model['ampl_season'][model_name]

#             # use precip weighted mean 
#             ampl_weighted_mean = ampl_weighted_mean_df[model_name]
#             phase_mean = circmean_phase_by_model[model_name]
#         else:
#             ecs_i = np.nan
#             ampl_mean = np.nan
#             phase_mean = np.nan
    
#     model_error_stats[model_name] = [model_i_std, model_i_corr, rmse_i, 
#                                      ampl_mean, ampl_weighted_mean, phase_mean, ecs_i]
    
# model_error_stats_df = pd.DataFrame(model_error_stats).T
# model_error_stats_df.columns = ['std', 'corr', 'rmse', 
#                                 'ampl_mean', 'ampl_weighted_mean', 'phase_mean', 'ecs']
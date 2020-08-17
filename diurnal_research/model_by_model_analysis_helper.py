'''Functions for computing stats on diurnal cycle indicies after being computed'''

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import circmean, circstd, mode

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

# def mode_apply_phase(df, round_precision = 1):
#     '''Handle binning of data before taking mode for given variables'''
#         df.round(round_precision)
# def mode_apply_ampl(df, round_precision = 4):


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
def _get_mean_field(data_dir,
                    field_name = 'rlut',
                    file_name = '1985-01_2006-01_mean.nc',
                    filter_lat = False,
                    landsea_mask_df = None,
                    landsea_bool = None,
                    **kwargs):
    '''Given path to cmip files (with subdirs for each model),
    take mean of each model and return dataframe.
    Args
    -----
    data_dir - str
        dir containing cmip models in folders
        
    filter_lat - bool
        whether or not to filter by lat. 
    landsea_mask_df - pd.DataFrame
        Pandas dataframe containing lat/lon as multiindex keys and 
        mask as column (0 for water, 1 for land). df should only
        contain one column. ``
    landsea_bool - bool 
        If 1 return land only, if 0 return sea only. 
            
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
                
                df_field = df_field_mask.loc[:,field_name]
                
            if filter_lat:
                df_field_filtered = filter_by_lat(df_field, **kwargs)
                var_dict[model_name] = df_field_filtered[field_name].mean()
#                 return df_field_filtered
            else:
                var_dict[model_name] = ds_field[field_name].mean().item()
            

        except Exception as e: 
            print('Could not process ', model_name, ' ', str(e))
    out_df = pd.DataFrame.from_dict(var_dict, orient = 'index')
    out_df.columns = [field_name,] 
    return out_df

def _get_mean_fields(
    field_names,
    mean_fields_to_rel_path_map,
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


# def fetch_means_n_35:
    
    
    


# def fetch_means_s_35:
    
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
    and return dataframe. Circular statistics are used for the phase field.
    
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
    
    _variables_of_interest = ['ampl_season', 'phase_season', 'model_name']
    
    if 'season' in df_for_stats.columns:
        _variables_of_interest.append('season')
    
    phase_bin_precision = 1
    ampl_bin_precision = 4
    
    if agg_method == 'mean':
    #     mean_field_by_model = df_for_stats.groupby('model_name').mean()
        agg_field_by_model = df_for_stats.groupby(groupby_vars).mean()
    #     circmean_phase_by_model = df_for_stats[['phase_season', 'model_name']].groupby('model_name').apply(phase_circmean)
        agg_phase_by_model = df_for_stats[_variables_of_interest].groupby(groupby_vars).agg(phase_circmean)
    
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
        else:
            rmse_i = np.nan
            model_i_corr = np.nan


        

        

        model_i_std = circstd(df_i[field].values, low = 0.0, high = 24.0)


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
    
    field_id_to_xlabel = {'phase_season': 'Phase [hours]',
                          'ampl_season': 'Amplitude'}
    if 'cmip_identifier' in df.columns:
        sns.distplot(df[(df['land_sea_mask'] == 0) & (df['cmip_identifier'] == cmip_identifier)][field_id].values, label = cmip_identifier + ' ' + 'Water', bins = nbins, ax = ax)
        sns.distplot(df[(df['land_sea_mask'] == 1) & (df['cmip_identifier'] == cmip_identifier)][field_id].values, label = cmip_identifier + ' ' + 'Land', bins = nbins, ax = ax)
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



def plot_corr_matrix(corr_mat_ds):
    plt.figure(figsize = (12,10))
    plt.title('Correlation Matrix for All Latitudes [Method: mode]')
    upper_tr_mask = np.triu(corr_mat_ds.corr())
    sns.heatmap(corr_mat_ds.corr(), annot = True, 
                vmin = -1, vmax = 1, center = 0, fmt='.2g',
                mask = upper_tr_mask)
    plt.tight_layout()
    
    

##### More plotting 


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
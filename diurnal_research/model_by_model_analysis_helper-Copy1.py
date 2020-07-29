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

def sin_hour(hours):
    '''Give array of hours, convert sin quanitity to deal with cirular quantity'''
    return np.sin(hours*HOURS_TO_RADIANS)

def ampl_weighted_mean_func(df):
    return (df['precip_weights']*df['ampl_season']).sum()

def mode_apply(df):
    '''Scipy mode returns both mode and counts, os wrap function to only return mode. '''
    return mode(df)[0].item()


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
    
    df = df.reset_index(['lat','lon'])
    
    if absolute_value:
        return df[(abs(df['lat'])>= min_lat) & (abs(df['lat']) <= max_lat)]

    else:
        return df[(df['lat']>= min_lat) & (df['lat'] <= max_lat)]

def filter_by_season(df, season):
    '''Filter df by season. '''
    df = df.reset_index('season')
    return df[df['season'] == season]
    
    
####### ds operators 
def _get_mean_field(data_dir,
                    field_name = 'rlut',
                    file_name = '1985-01_2006-01_mean.nc',
                    filter_lat = False,
                    **kwargs):
    '''Given path to cmip files (with subdirs for each model),
    take mean of each model and return dict.
    Args
    -----
    data_dir - str
        dir containing cmip models in folders
        
        filter_lat - bool
            whether or not to filter by lat. 
            
        **kwargs
            keyword arguments to be past to lat filtering function
            (filter_by_lat)
    
    '''
    rad_dict = {}

    model_names = os.listdir(data_dir)
    for model_name in model_names:
        try:
            ds_field = xr.open_dataset(data_dir + model_name + '/' + file_name)
            
            if filter_lat:
                df_field = ds_field.to_dataframe()
                df_field_filtered = filter_by_lat(df_field, **kwargs)
                rad_dict[model_name] = df_field_filtered[field_name].mean()
                
            else:
                rad_dict[model_name] = ds_field[field_name].mean().item()
            

        except Exception as e: 
            print('Could not process ', model_name, ' ', str(e))
    
    return rad_dict

    
    
# compute stats on df containing indicies 

def compute_stats(df_for_stats,
                 df_for_stats_true,
                 field = 'phase_season',
                 agg_method = 'mean',
                 additional_stats = True, 
                 ecs_dict = None,
                 tcr_dict = None,
                 rlut_dict = None,
                 rsut_dict = None):
    '''Given pd.DataFrames containing cmip model and satellite diurnal cycle indicies, compute various error statistics
    and return dataframe. Circular statistics are used for the phase field.
    
    Args
    -----------
        df_for_stats - pd.DataFrame
            Dataframe containing CMIP model diurnal stats
        df_for_stats_true - pd.DataFrame
            Dataframe containing data to validate agaisnt. 
        
        field - str
        
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
    
    # subset columns before doing pandas aggreations to speed up calculations
    _variables_of_interest = ['ampl_season', 'model_name']
    phase_bin_precision = 1
    ampl_bin_precision = 3
    
    if agg_method == 'mean':
    #     mean_field_by_model = df_for_stats.groupby('model_name').mean()
        agg_field_by_model = df_for_stats.groupby('model_name').mean()
    #     circmean_phase_by_model = df_for_stats[['phase_season', 'model_name']].groupby('model_name').apply(phase_circmean)
        agg_phase_by_model = df_for_stats[['phase_season', 'model_name']].groupby('model_name').agg(phase_circmean)
    
    elif agg_method == 'mode':
        agg_field_by_model = df_for_stats[_variables_of_interest].round(ampl_bin_precision).groupby('model_name').agg(mode_apply)
        
        agg_phase_by_model = df_for_stats[['phase_season', 'model_name']].round(phase_bin_precision).groupby('model_name').agg(mode_apply)

    for model_name, df_i in df_for_stats.groupby('model_name'):
        
        rmse_i = circrmse(df_i[field], df_for_stats_true[field])

        model_i_corr = np.corrcoef(sin_hour(df_i[field].values), sin_hour(df_for_stats_true[field].values))[0,1]
    #     model_i_corr = np.corrcoef(df_i[field].values, df_gpm[field].values)[0,1]
    #     model_i_corr = astro_circstats.circcorrcoef(df_i[field].values * HOURS_TO_RADIANS,  
    #                              df_gpm[field].values * HOURS_TO_RADIANS)

        model_i_std = circstd(df_i[field].values, low = 0.0, high = 24.0)


        # To Do
        if additional_stats == True:
         # use regular mean
#             ampl_mean = mean_field_by_model['ampl_season'][model_name]
            ampl_mean = agg_field_by_model['ampl_season'][model_name]

            # use precip weighted mean 
#                 ampl_weighted_mean = ampl_weighted_mean_df[model_name]
            phase_mean = agg_phase_by_model['phase_season'][model_name]
            
        
            if (not ecs_dict is None) & (model_name in ecs_dict):
                ecs_i = ecs_dict[model_name]
            else:
                ecs_i = np.nan
                
            if (not tcr_dict is None) & (model_name in tcr_dict):
                tcr_i = tcr_dict[model_name]
            else:
                tcr_i = np.nan
                
            if (not rlut_dict is None) & (model_name in rlut_dict):
                rlut_i = rlut_dict[model_name]
            else:
                rlut_i = np.nan
            
            if (not rsut_dict is None) & (model_name in rsut_dict):
                rsut_i = rsut_dict[model_name]
            else:
                rsut_i = np.nan
                

        model_error_stats[model_name] = [model_i_std, model_i_corr, rmse_i, 
                                         ampl_mean, phase_mean, ecs_i, tcr_i,
                                         rlut_i, rsut_i]

    model_error_stats_df = pd.DataFrame(model_error_stats).T
    model_error_stats_df.columns = ['std', 'corr', 'rmse', 
                                    'ampl_mean', 'phase_mean', 'ecs', 'tcr',
                                    'rlut', 'rsut']
    return model_error_stats_df

################ Plotting #############################
def land_sea_histogram(df, 
                       nbins = 20,
                       cmip_identifier = 'CMIP6', 
                       field_id = 'phase_season'):
    
    field_id_to_xlabel = {'phase_season': 'Phase [hours]',}
    
    sns.distplot(df[(df['land_sea_mask'] == 0) & (df['cmip_indentifier'] == cmip_identifier)][field_id].values, label = 'Water', bins = nbins)
    sns.distplot(df[(df['land_sea_mask'] == 1) & (df['cmip_indentifier'] == cmip_identifier)][field_id].values, label = 'Land', bins = nbins)

    plt.xlim([0, 24])
    plt.ylim([0, 0.25])
    plt.xlabel(field_id_to_xlabel[field_id])
    plt.ylabel('Probability Density')
    plt.title('PDF of CMIP6 Diurnal Precipitation Phase [hours]')
    plt.grid()
    plt.legend()
    
    return plt.gca()

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
'''Functions for computing stats on diurnal cycle indicies after being computed'''

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import circmean, circstd

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



# compute stats on df containing indicies 

def compute_stats(df_for_stats,
                 df_for_stats_true,
                 field = 'phase_season',
                 additional_stats = True, 
                 ecs_dict = None,
                 tcr_dict = None):
    '''Given pd.DataFrames containing cmip model and satellite diurnal cycle indicies, compute various error statistics
    and return dataframe. Circular statistics are used for the phase field.
    
    Args
    -----------
        df_for_stats - pd.DataFrame
            Dataframe containing CMIP model diurnal stats
        df_for_stats_true - pd.DataFrame
            Dataframe containing data to validate agaisnt. 
        
        field - str
        
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
    

    mean_field_by_model = df_for_stats.groupby('model_name').mean()
    circmean_phase_by_model = df_for_stats[['phase_season', 'model_name']].groupby('model_name').apply(phase_circmean)


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
            ampl_mean = mean_field_by_model['ampl_season'][model_name]

            # use precip weighted mean 
#                 ampl_weighted_mean = ampl_weighted_mean_df[model_name]
            phase_mean = circmean_phase_by_model[model_name]
    
            if model_name in ecs_dict:
                ecs_i = ecs_dict[model_name]
            else:
                ecs_i = np.nan
                
            if model_name in tcr_dict:
                tcr_i = tcr_dict[model_name]
            else:
                tcr_i = np.nan
#                 ampl_mean = np.nan
#                 phase_mean = np.nan

        model_error_stats[model_name] = [model_i_std, model_i_corr, rmse_i, 
                                         ampl_mean, phase_mean, ecs_i, tcr_i]

    model_error_stats_df = pd.DataFrame(model_error_stats).T
    model_error_stats_df.columns = ['std', 'corr', 'rmse', 
                                    'ampl_mean', 'phase_mean', 'ecs', 'tcr']
    return model_error_stats_df



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
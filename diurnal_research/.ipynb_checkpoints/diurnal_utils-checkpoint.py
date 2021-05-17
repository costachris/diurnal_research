import os
import cftime
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ephem

import pickle
import ipdb

import math
from datetime import datetime, time, timedelta
from math import pi, cos, sin
from scipy import optimize
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
import numba
from joblib import Parallel, delayed


HOURS_TO_RADIANS = 2*np.pi/24
### Utility functions for finding amp, phase of diurnal cycle

def cos_func(x, a, phi):
    b = ((2*np.pi)/24)
    return a * np.cos(b * (x - phi))

def cos_func_2_mode(x, a1, phi1, a2, phi2):
    b = ((2*np.pi)/24)
    c = ((2*np.pi)/12)
    return (a1 * np.cos(b * (x - phi1)) + a2 * np.cos(c * (x - phi2)))

def cos_fit_grid(field_array, lst_array, hour_means, hour_bins):
    ''' Given ds, perform cos fit grid-by-grid and return array of fitted parameters.
        Inputs: 
            field_array - [time, lat, lon]
            lst_array = [time, lat, lon]
    
        Returns: 
            amplitude_xy - [lat, lon]
            frequency_xy - [lat, lon]
            phase_xy - [lat, lon] 
            '''
    
    out_shape = (field_array.shape[1], field_array.shape[2])
    amplitude_xy = np.zeros(out_shape)
    frequency_xy = np.zeros(out_shape)
    phase_xy = np.zeros(out_shape)
    
    for lat_ii in tqdm(range(out_shape[0])):
        for lon_ii in range(out_shape[1]):
            
            ts_loc = field_array[:,lat_ii, lon_ii]
            time_loc = lst_array[:,lat_ii, lon_ii]
            ts_loc_hourly_mean = hour_means[:,lat_ii, lon_ii]
            
            hours_to_radians = 2*np.pi/24
            try:

                params, params_covariance = optimize.curve_fit(cos_func, 
                                                               time_loc, 
                                                               ts_loc - ts_loc.mean(),
                                                               bounds = (0,24),
                                                               p0=[10, 
                                                                   hour_bins[np.argmax(ts_loc_hourly_mean)]],
                                                               maxfev=10000)
            except: 
                params = np.array([np.nan, np.nan, np.nan])

            amplitude_xy[lat_ii, lon_ii] = params[0]
            phase_xy[lat_ii, lon_ii] = params[1]
        
    return(amplitude_xy, phase_xy)


def cos_fit_2_mode_grid_average(field_array, hour_bins):
    ''' Given array of timeseries along time axis, perform cos fit 
    grid-by-grid and return array of fitted phase and ampltiude. Returns
    grid of fitted parameters and associated uncertainties.
    
        Inputs: 
            field_array - [time, lat, lon]
                Array of values to perform cos fit for. 
            hour_bins - [1xtime] array
                
                
        Returns: 
            amplitude_xy - [lat, lon]
                Amplitude of fit. 
                
            phase_xy - [lat, lon]
                Phase of fit.
            amplitude_xy_cov - 
                Covariance of amplitude fit
            phase_xy-cov - 
                Covariance of phase fit
            '''
    
    out_shape = (field_array.shape[1], field_array.shape[2])
    amplitude_xy = np.zeros(out_shape)
    phase_xy = np.zeros(out_shape)
    
    # cov matricies
    amplitude_xy_cov = np.zeros(out_shape)
    phase_xy_cov = np.zeros(out_shape)
    
#     tt1 = np.zeros(out_shape)
#     tt2 = np.zeros(out_shape)
    
    for lat_ii in tqdm(range(out_shape[0])):
        for lon_ii in range(out_shape[1]):
            
            ts_loc = field_array[:,lat_ii, lon_ii]
#             ts_loc_hourly_mean = field_array[:,lat_ii, lon_ii]
#             ipdb.set_trace()
            try:

                ts_loc_mean_sub = ts_loc - np.nanmean(ts_loc)
                non_nan_inds = np.where(np.isfinite(ts_loc_mean_sub))

                ts_loc_mean_sub_interp = ts_loc[non_nan_inds]
                hour_bins_non_nan =  hour_bins[non_nan_inds]

                params, params_covariance = optimize.curve_fit(cos_func_2_mode, 
                                                               hour_bins_non_nan, 
                                                               ts_loc_mean_sub_interp,
#                                                                p0=[ts_loc_mean_sub_interp.std(),
#                                                                    hour_bins_non_nan[np.nanargmax(ts_loc_mean_sub_interp)]],
                                                               p0=[ts_loc_mean_sub_interp.std(),
                                                                   hour_bins_non_nan[np.nanargmax(ts_loc_mean_sub_interp)],
                                                                   ts_loc_mean_sub_interp.std(),
                                                                   hour_bins_non_nan[np.nanargmax(ts_loc_mean_sub_interp)],
                                                                  
                                                                  ],
                                                                   maxfev=10000000)
            except: 
                params = np.array([np.nan, np.nan])
                params_covariance = np.array([[np.nan, np.nan], [np.nan, np.nan]])
        
            

            amplitude_xy[lat_ii, lon_ii] = params[0]
            
            # the %24 below acccounts for cases when phase is outside of [0, 24] hours
            phase_xy[lat_ii, lon_ii] = params[1]
            
#             tt1[lat_ii, lon_ii], tt2[lat_ii, lon_ii] = params[2], params[3]
            # handle cases where amplitude is negative by making postive and
            # shifting phase by pi. 
            if amplitude_xy[lat_ii, lon_ii] < 0:
                amplitude_xy[lat_ii, lon_ii] = -1*amplitude_xy[lat_ii, lon_ii]
                phase_xy[lat_ii, lon_ii] += HOURS_TO_RADIANS*12
            
#             return params_covariance
            phase_xy[lat_ii, lon_ii] = phase_xy[lat_ii, lon_ii]%24
            # cov measures
            amplitude_xy_cov[lat_ii, lon_ii] = np.diag(params_covariance)[0]
            phase_xy_cov[lat_ii, lon_ii] = np.diag(params_covariance)[1]
        
    return(amplitude_xy, phase_xy, amplitude_xy_cov, phase_xy_cov)




def cos_fit_grid_average(field_array, hour_bins):
    ''' Given array of timeseries along time axis, perform cos fit 
    grid-by-grid and return array of fitted phase and ampltiude. Returns
    grid of fitted parameters and associated uncertainties.
    
        Inputs: 
            field_array - [time, lat, lon]
                Array of values to perform cos fit for. 
            hour_bins - [1xtime] array
                
                
        Returns: 
            amplitude_xy - [lat, lon]
                Amplitude of fit. 
                
            phase_xy - [lat, lon]
                Phase of fit.
            amplitude_xy_cov - 
                Covariance of amplitude fit
            phase_xy-cov - 
                Covariance of phase fit
            '''
    
    out_shape = (field_array.shape[1], field_array.shape[2])
    amplitude_xy = np.zeros(out_shape)
    phase_xy = np.zeros(out_shape)
    
    # cov matricies
    amplitude_xy_cov = np.zeros(out_shape)
    phase_xy_cov = np.zeros(out_shape)
    
    for lat_ii in tqdm(range(out_shape[0])):
        for lon_ii in range(out_shape[1]):
            
            ts_loc = field_array[:,lat_ii, lon_ii]
#             ts_loc_hourly_mean = field_array[:,lat_ii, lon_ii]
            try:

                ts_loc_mean_sub = ts_loc - np.nanmean(ts_loc)
                non_nan_inds = np.where(np.isfinite(ts_loc_mean_sub))

                f = interpolate.interp1d(hour_bins[non_nan_inds], 
                                         ts_loc_mean_sub[non_nan_inds], 
                                         kind = 'quadratic',
                                         fill_value="extrapolate")

                ts_loc_mean_sub_interp = f(hour_bins)

                params, params_covariance = optimize.curve_fit(cos_func, 
                                                               hour_bins, 
                                                               ts_loc_mean_sub_interp,
#                                                                bounds = (0,24),
                                                               p0=[ts_loc_mean_sub_interp.std(),
                                                                   hour_bins[np.nanargmax(ts_loc)]],
                                                               maxfev=10000)
            except: 
                params = np.array([np.nan, np.nan])

            amplitude_xy[lat_ii, lon_ii] = params[0]
            
            # the %24 below acccounts for cases when phase is outside of [0, 24] hours
            phase_xy[lat_ii, lon_ii] = params[1]
            
            # handle cases where amplitude is negative by making postive and
            # shifting phase by pi. 
            if amplitude_xy[lat_ii, lon_ii] < 0:
                amplitude_xy[lat_ii, lon_ii] = -1*amplitude_xy[lat_ii, lon_ii]
                phase_xy[lat_ii, lon_ii] += HOURS_TO_RADIANS*12
            
            phase_xy[lat_ii, lon_ii] = phase_xy[lat_ii, lon_ii]%24
            # cov measures
            amplitude_xy_cov[lat_ii, lon_ii] = np.diag(params_covariance)[0]
            phase_xy_cov[lat_ii, lon_ii] = np.diag(params_covariance)[1]
        
    return(amplitude_xy, phase_xy, amplitude_xy_cov, phase_xy_cov)



### Utility function for calculating local solar time 
def calc_solar_time(utc_time, 
                    longitude, 
                    bin_interval = 3, 
                    round_precision = 2, 
                    bin_bool = True):
    ''' Compute local solartime give longitude and time array-like objects.
    Args - 
        bin_interval - int 
    
    '''
    longitude_radians = np.radians(longitude)
#     return utc_time + timedelta(hours=longitude_radians / np.pi * 12)
    lst_approx = ((utc_time.hour + utc_time.minute/60.0) + (longitude_radians / np.pi * 12))#%24
    if bin_bool:
#         return (lst_approx - (lst_approx%3))
        binned_lst = ((lst_approx + (bin_interval/2)) - ((lst_approx+(bin_interval/2))%bin_interval))%24
        return np.round(binned_lst, round_precision)
    else:
        return lst_approx%24

def calc_solar_time_ephem_single_grid(dt, 
                          lon, 
                          bin_interval = 3, 
                          bin_bool = True):
    '''
    Args: 
        lat - lattitude in degrees
        lon -longitude in degrees
    
    '''
    obs = ephem.Observer()
#     obs.date = dt.strftime('%Y/%m/%d %H:%M:%S')
    obs.date = dt
    obs.lon = np.radians(lon)
#     obs.lat = np.radians(lat)
#     obs.elevation = elev
    sun = ephem.Sun()
    
    sun.compute(obs)
    hour_angle = obs.sidereal_time() - sun.ra
#     return hour_angle
    lst_approx = (ephem.hours(hour_angle + ephem.hours('12:00')).norm)*(24/(2*np.pi))  # norm for 24h
    
    if bin_bool:
        return ((lst_approx + (bin_interval/2)) - ((lst_approx+(bin_interval/2))%bin_interval))%24
    else:
        return lst_approx%24

# @numba.jit()
def calc_solar_time_ephem(dt, 
                          lon, 
                          bin_interval = 3, 
                          parallel = False, 
                          bin_bool = True):
    '''Uses calc_solar_time_ephem_single_grid function to compute LST over grid of lons, given UTC time.'''
    lst_approx = np.zeros(lon.shape)
    
    if parallel:
        lst_approx = Parallel(n_jobs=40, verbose = 0)(delayed(calc_solar_time_ephem_single_grid)(dt, 
                                                              lon[lat_ii, lon_ii],
                                                              bin_interval = bin_interval,
                                                              bin_bool = bin_bool) for lat_ii in range(lon.shape[0]) for lon_ii in range(lon.shape[1]))
   
    else:
        for lat_ii in range(lon.shape[0]):
            for lon_ii in range(lon.shape[1]):
                lst_approx[lat_ii, lon_ii] = calc_solar_time_ephem_single_grid(dt, 
                                                              lon[lat_ii, lon_ii],
                                                              bin_interval = bin_interval,
                                                              bin_bool = bin_bool)
    return lst_approx

    
def compute_lst_array(ds, 
                      bin_interval = 3, 
                      bin_bool = True, 
#                       method = 'ephem',
                      method = 'simple',
                      field_id = 'clt', 
                      lon_id = 'lon', 
                      lat_id = 'lat', 
                      time_id = 'time', 
                      round_precision = 2, 
                      lon_mesh = None, 
                      lat_mesh = None):
    '''Given xr.dataset with longitude and UTC time field, compute lst array (time, lat, lon).
        Relies on calc_solar_time function. '''
    if (lon_mesh is None) and (lat_mesh is None):
        lon_mesh, lat_mesh = np.meshgrid(ds[lon_id].values, ds[lat_id].values)
    
    lst_times = np.empty(ds[field_id].shape).astype(np.float32)

    for time_ind in tqdm(range(len(ds[time_id]))):
        # note: for diurnal analysis we don't really care about year, month, day: just hour
        if (type(ds['time'][0].item()) == cftime._cftime.DatetimeNoLeap) | \
            (type(ds['time'][0].item()) == cftime._cftime.Datetime360Day) | \
            (type(ds['time'][0].item()) == cftime._cftime.DatetimeJulian):
#             dt_i = datetime.strptime(str(ds.isel(time = time_ind)[time_id].values.item()), '%Y-%m-%d %H:%M:%S')
            dt_i = ds.isel(time = time_ind)[time_id].values.item()
            
        else:
            dt_i = datetime.utcfromtimestamp(ds[time_id][time_ind].item() * 1e-9)
        
        if method == 'simple':
            lst_times[time_ind,:,:] = calc_solar_time(dt_i, 
                                                      lon_mesh,
                                                      bin_interval = bin_interval,
                                                      round_precision = round_precision, 
                                                      bin_bool = bin_bool)
        if method == 'ephem': 
#             print('Using ephem method to calculate LST.')
            lst_times[time_ind,:,:] = calc_solar_time_ephem(dt_i, 
                                                            lon_mesh, 
                                                            bin_interval = bin_interval, 
                                                            bin_bool = bin_bool)
#             print('here', lst_times.dtype)
    lst_da = xr.DataArray(lst_times, dims = (time_id, lat_id, lon_id), coords = {time_id:ds[time_id].values,
                                                                        lon_id: ds[lon_id].values, 
                                                                        lat_id: ds[lat_id].values})
    return lst_da


#########
def diurnal_analysis_2_mode_full_bootstrapped_int_days_precomputed(
                     ds,
                     lst_array_full, 
                     field_season_array_full,
                     grid_time_resolution_hours = 3, 
                     time_resolution_hours = 1,
                     round_precision = 2,
                     **bootstrap_kwargs):
    ''' Perform diurnal analysis using stationary bootstrapping (integer number of days) to get uncertainty on parameters. A lot of arrays are precomputed (LST, field_array) to speed up 
    bootstrap calculations. 
    
    Assumes ordering (time, lat, lon), with those field/coordinate names.
    Given dataset compute mean, standard deviation, amplitude, and phase. '''
    hour_bins = np.arange(time_resolution_hours, 24 + time_resolution_hours, time_resolution_hours)
    hour_bins = np.round(hour_bins, round_precision)
    grid_hour_bins = np.arange(grid_time_resolution_hours, 24 + grid_time_resolution_hours, grid_time_resolution_hours)
    grid_hour_bins = np.round(grid_hour_bins, round_precision)
    lon_mesh, lat_mesh = np.meshgrid(ds['lon'].values, ds['lat'].values)

    mu_season = {}
    sigma_season = {}
    ampl_season = {}
    phase_season = {}
    average_cycle_season = {}
    
    datasize = len(lst_array_full)
    bootstrap_inds = stationary_bootstrap_int(datasize, **bootstrap_kwargs)
    print('bootstrap_shape', bootstrap_inds.shape)
    for bootstrap_number in range(bootstrap_inds.shape[0]):
        
        print(bootstrap_number)
        # choose subset of dataset using bootstrap indicies
        bootstrap_inds_i = bootstrap_inds[bootstrap_number]
        field_season_array =  field_season_array_full[bootstrap_inds_i,:,:]
        lst_array = lst_array_full[bootstrap_inds_i,:,:]
        # compute hourly grid means needed for cos fit
        f_bar_ks = {}
        for ii in range(len(hour_bins)):
            print(ii)
            hour_i = hour_bins[ii]
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_ks[hour_i] = f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))

        print('Performing Cos Fit')
        res = cos_fit_2_mode_grid_average(hour_means, hour_bins)
        print('Finished Cos Fit')
        ampl_season[bootstrap_number], phase_season[bootstrap_number] = res[0], res[1]
        
    out_ds = xr.Dataset()
#     out_ds['mu_season'] = make_da_from_dict(mu_season, ds, time_id = 'bootstrap_id')
    out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds, time_id = 'bootstrap_id')
    out_ds['phase_season'] = make_da_from_dict(phase_season,ds, time_id = 'bootstrap_id')

#     out_ds_means = xr.Dataset()
#     out_ds_means = make_da_from_dict_time(average_cycle_season, ds, grid_hour_bins)
#     out_ds_means = out_ds_means.to_dataset(name = field_id + '_mean')
#     return (out_ds, out_ds_means)
    return out_ds

def diurnal_analysis_2_mode_full_bootstrapped_int_days(ds, 
                     field_id, 
                     grid_time_resolution_hours = 3, 
                     time_resolution_hours = 1,
                     round_precision = 2,
                     **bootstrap_kwargs):
    ''' Perform diurnal analysis using stationary bootstrapping (integer number of days) to get uncertainty on parameters. 
    Assumes ordering (time, lat, lon), with those field/coordinate names.
    Given dataset compute mean, standard deviation, amplitude, and phase. '''
    hour_bins = np.arange(time_resolution_hours, 24 + time_resolution_hours, time_resolution_hours)
    hour_bins = np.round(hour_bins, round_precision)
    grid_hour_bins = np.arange(grid_time_resolution_hours, 24 + grid_time_resolution_hours, grid_time_resolution_hours)
    grid_hour_bins = np.round(grid_hour_bins, round_precision)
    lon_mesh, lat_mesh = np.meshgrid(ds['lon'].values, ds['lat'].values)

    mu_season = {}
    sigma_season = {}
    ampl_season = {}
    phase_season = {}
    average_cycle_season = {}
    
    datasize = len(ds['time'])
    bootstrap_inds = stationary_bootstrap_int(datasize, **bootstrap_kwargs)
    print('bootstrap_shape', bootstrap_inds.shape)
    for bootstrap_number in range(bootstrap_inds.shape[0]):
        
        print(bootstrap_number)
        # choose subset of dataset using bootstrap indicies
        ds_sub = ds.isel(time = bootstrap_inds[bootstrap_number])
#         lst_da = compute_lst_array(ds_sub, 
#                                    bin_interval = grid_time_resolution_hours, # time_resolution_hours,
#                                    bin_bool = True, 
#                                    round_precision = round_precision,
#                                    lon_mesh = lon_mesh, 
#                                    lat_mesh = lat_mesh,
#                                    field_id = field_id)
#         lst_array = lst_da.values
#         lst_array = lst_array.astype(np.float32)
#         del lst_da
#         lst_array[lst_array == 0] = 24

        field_season_array = ds_sub[field_id].values
        field_season_array =  field_season_array.astype(np.float32)


#         mu_ij = np.zeros(field_season_array.shape[-2:])

#         f_bar_ks = {}
#         for hour_i in grid_hour_bins:
#             masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

#             # mean for a given season, LST
#             f_bar_k = np.nanmean(masked_field, axis = 0)
#             f_bar_k[np.isnan(f_bar_k)] = 0
#             f_bar_ks[hour_i] = f_bar_k
#             mu_ij  += f_bar_k

#         hour_means = np.stack(list(f_bar_ks.values()))
#         average_cycle_season[bootstrap_number] = hour_means


#         mu_ij = (1/len(grid_hour_bins))*mu_ij
#         mu_season[bootstrap_number] = mu_ij


        lst_da = compute_lst_array(ds_sub, 
                                   bin_interval = time_resolution_hours,
                                   bin_bool = True, 
                                   round_precision = round_precision,
                                   lon_mesh = lon_mesh, 
                                   lat_mesh = lat_mesh,
                                   field_id = field_id)
        lst_array = lst_da.values
        lst_array = lst_array.astype(np.float32)
        del lst_da
        lst_array[lst_array == 0] = 24
        field_season_array = ds_sub[field_id].values
        field_season_array = field_season_array.astype(np.float32)
        
        # compute hourly grid means needed for cos fit
        f_bar_ks = {}
        for ii in range(len(hour_bins)):
            hour_i = hour_bins[ii]
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_ks[hour_i] = f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))

        print('Performing Cos Fit')
        res = cos_fit_2_mode_grid_average(hour_means, hour_bins)
        print('Finished Cos Fit')
        ampl_season[bootstrap_number], phase_season[bootstrap_number] = res[0], res[1]
        
    out_ds = xr.Dataset()
#     out_ds['mu_season'] = make_da_from_dict(mu_season, ds, time_id = 'bootstrap_id')
    out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds, time_id = 'bootstrap_id')
    out_ds['phase_season'] = make_da_from_dict(phase_season,ds, time_id = 'bootstrap_id')

#     out_ds_means = xr.Dataset()
#     out_ds_means = make_da_from_dict_time(average_cycle_season, ds, grid_hour_bins)
#     out_ds_means = out_ds_means.to_dataset(name = field_id + '_mean')
#     return (out_ds, out_ds_means)
    return out_ds



def diurnal_analysis_2_mode_full_bootstrapped(ds, 
                     field_id, 
                     grid_time_resolution_hours = 3, 
                     time_resolution_hours = 1,
                     round_precision = 2,
                     **bootstrap_kwargs):
    ''' Perform diurnal analysis using stationary bootstrapping to get uncertainty on parameters. 
    Assumes ordering (time, lat, lon), with those field/coordinate names.
    Given dataset compute mean, standard deviation, amplitude, and phase. '''
    hour_bins = np.arange(time_resolution_hours, 24 + time_resolution_hours, time_resolution_hours)
    hour_bins = np.round(hour_bins, round_precision)
    grid_hour_bins = np.arange(grid_time_resolution_hours, 24 + grid_time_resolution_hours, grid_time_resolution_hours)
    grid_hour_bins = np.round(grid_hour_bins, round_precision)
    lon_mesh, lat_mesh = np.meshgrid(ds['lon'].values, ds['lat'].values)

    mu_season = {}
    sigma_season = {}
    ampl_season = {}
    phase_season = {}
    average_cycle_season = {}
    
    datasize = len(ds['time'])
    bootstrap_inds = stationary_bootstrap(datasize, **bootstrap_kwargs)
#     return bootstrap_inds
    for bootstrap_number in range(bootstrap_inds.shape[0]):
        
        print(bootstrap_number)
        # choose subset of dataset using bootstrap indicies
        ds_sub = ds.isel(time = bootstrap_inds[bootstrap_number])
        lst_da = compute_lst_array(ds_sub, 
                                   bin_interval = grid_time_resolution_hours, # time_resolution_hours,
                                   bin_bool = True, 
                                   round_precision = round_precision,
                                   lon_mesh = lon_mesh, 
                                   lat_mesh = lat_mesh,
                                   field_id = field_id)
        lst_array = lst_da.values
        lst_array = lst_array.astype(np.float32)
        del lst_da
        lst_array[lst_array == 0] = 24

        field_season_array = ds_sub[field_id].values
        field_season_array =  field_season_array.astype(np.float32)


        mu_ij = np.zeros(field_season_array.shape[-2:])

        f_bar_ks = {}
        for hour_i in grid_hour_bins:
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_k[np.isnan(f_bar_k)] = 0
            f_bar_ks[hour_i] = f_bar_k
            mu_ij  += f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))
        average_cycle_season[bootstrap_number] = hour_means


        mu_ij = (1/len(grid_hour_bins))*mu_ij
        mu_season[bootstrap_number] = mu_ij


        lst_da = compute_lst_array(ds_sub, 
                                   bin_interval = time_resolution_hours,
                                   bin_bool = True, 
                                   round_precision = round_precision,
                                   lon_mesh = lon_mesh, 
                                   lat_mesh = lat_mesh,
                                   field_id = field_id)
        lst_array = lst_da.values
        del lst_da
        lst_array[lst_array == 0] = 24
        field_season_array = ds_sub[field_id].values

        # compute hourly grid means needed for cos fit
        f_bar_ks = {}
        for ii in range(len(hour_bins)):
            hour_i = hour_bins[ii]
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_ks[hour_i] = f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))

        print('Performing Cos Fit')
        res = cos_fit_2_mode_grid_average(hour_means, hour_bins)
        print('Finished Cos Fit')
        ampl_season[bootstrap_number], phase_season[bootstrap_number] = res[0], res[1]
        
    out_ds = xr.Dataset()
    out_ds['mu_season'] = make_da_from_dict(mu_season, ds, time_id = 'bootstrap_id')
    out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds, time_id = 'bootstrap_id')
    out_ds['phase_season'] = make_da_from_dict(phase_season,ds, time_id = 'bootstrap_id')

    out_ds_means = xr.Dataset()
    out_ds_means = make_da_from_dict_time(average_cycle_season, ds, grid_hour_bins)
    out_ds_means = out_ds_means.to_dataset(name = field_id + '_mean')
       

    return (out_ds, out_ds_means)



def diurnal_analysis_2_mode_full(ds, 
                     field_id, 
                     grid_time_resolution_hours = 3, 
                     time_resolution_hours = 1,
                     round_precision = 2):
    ''' Assumes ordering (time, lat, lon), with those field/coordinate names.
    Given dataset compute mean, standard deviation, amplitude, and phase. '''
    hour_bins = np.arange(time_resolution_hours, 24 + time_resolution_hours, time_resolution_hours)
    hour_bins = np.round(hour_bins, round_precision)
    grid_hour_bins = np.arange(grid_time_resolution_hours, 24 + grid_time_resolution_hours, grid_time_resolution_hours)
    grid_hour_bins = np.round(grid_hour_bins, round_precision)
    lon_mesh, lat_mesh = np.meshgrid(ds['lon'].values, ds['lat'].values)

    mu_season = {}
    sigma_season = {}
    ampl_season = {}
    phase_season = {}
    ampl_cov_season = {}
    phase_cov_season = {}
    average_cycle_season = {}

    season_i = 'year_mean'

    lst_da = compute_lst_array(ds, 
                               bin_interval = grid_time_resolution_hours, # time_resolution_hours,
                               bin_bool = True, 
                               round_precision = round_precision,
                               lon_mesh = lon_mesh, 
                               lat_mesh = lat_mesh,
                               field_id = field_id)
    lst_array = lst_da.values
    lst_array = lst_array.astype(np.float32)
    del lst_da
    lst_array[lst_array == 0] = 24

    field_season_array = ds[field_id].values
    field_season_array =  field_season_array.astype(np.float32)
#         field_season_mean = field_season_array.mean(axis = 0)
    print('float32')
    # compute mu_ij
#         mu_ij = np.zeros(field_season_array.shape[-2:])
    mu_ij = np.zeros(field_season_array.shape[-2:])
#         sigma_ij = np.zeros((len(hour_bins),) + field_season_array.shape[-2:])

    f_bar_ks = {}
    for hour_i in grid_hour_bins:
        masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

        # mean for a given season, LST
        f_bar_k = np.nanmean(masked_field, axis = 0)
        f_bar_k[np.isnan(f_bar_k)] = 0
        f_bar_ks[hour_i] = f_bar_k
        mu_ij  += f_bar_k

    hour_means = np.stack(list(f_bar_ks.values()))
    average_cycle_season[season_i] = hour_means


    mu_ij = (1/len(grid_hour_bins))*mu_ij
#         mu_ij = (1/8)*mu_ij
    mu_season[season_i] = mu_ij

    sigma_ij = np.zeros(field_season_array.shape[-2:])
    for hour_i in grid_hour_bins:
        masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

        # mean for a given season, LST
        f_bar_k = np.nanmean(masked_field, axis = 0)
        f_bar_k[np.isnan(f_bar_k)] = 0
        sigma_ij += np.square(f_bar_k - mu_ij)

    sigma_ij = np.sqrt((1/(len(grid_hour_bins) - 1))*sigma_ij)
#         sigma_ij = np.sqrt((1/7)*sigma_ij)
#         sigma_ij = np.nanstd(sigma_ij, axis = 0)

    sigma_season[season_i] = sigma_ij

    lst_da = compute_lst_array(ds, 
                               bin_interval = time_resolution_hours,
                               bin_bool = True, 
                               round_precision = round_precision,
                               lon_mesh = lon_mesh, 
                               lat_mesh = lat_mesh,
                               field_id = field_id)
    lst_array = lst_da.values
    del lst_da
    lst_array[lst_array == 0] = 24
    field_season_array = ds[field_id].values

    # compute hourly grid means needed for cos fit
    f_bar_ks = {}
    for ii in range(len(hour_bins)):
        hour_i = hour_bins[ii]
        masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

        # mean for a given season, LST
        f_bar_k = np.nanmean(masked_field, axis = 0)
        f_bar_ks[hour_i] = f_bar_k

    hour_means = np.stack(list(f_bar_ks.values()))

    print('Performing Cos Fit')


    res = cos_fit_2_mode_grid_average(hour_means, hour_bins)
    print('Finished Cos Fit')
    ampl_season[season_i], phase_season[season_i] = res[0], res[1]
    ampl_cov_season[season_i], phase_cov_season[season_i] = res[2], res[3]
        
    out_ds = xr.Dataset()
    out_ds['mu_season'] = make_da_from_dict(mu_season, ds)
    out_ds['sigma_season'] = make_da_from_dict(sigma_season,ds)
    out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds)
    out_ds['phase_season'] = make_da_from_dict(phase_season,ds)
    out_ds['ampl_cov_season'] = make_da_from_dict(ampl_cov_season, ds)
    out_ds['phase_cov_season'] = make_da_from_dict(phase_cov_season,ds)

    out_ds_means = xr.Dataset()
    out_ds_means = make_da_from_dict_time(average_cycle_season, ds, grid_hour_bins)
    out_ds_means = out_ds_means.to_dataset(name = field_id + '_mean')
       

    return (out_ds, out_ds_means)

def diurnal_analysis_2_mode_seasonal(ds, 
                     field_id, 
                     grid_time_resolution_hours = 3, 
                     time_resolution_hours = 1,
                     round_precision = 2):
    ''' Assumes ordering (time, lat, lon), with those field/coordinate names.
    Given dataset compute mean, standard deviation, amplitude, and phase. '''
    hour_bins = np.arange(time_resolution_hours, 24 + time_resolution_hours, time_resolution_hours)
    hour_bins = np.round(hour_bins, round_precision)
    grid_hour_bins = np.arange(grid_time_resolution_hours, 24 + grid_time_resolution_hours, grid_time_resolution_hours)
    grid_hour_bins = np.round(grid_hour_bins, round_precision)
    
    lon_mesh, lat_mesh = np.meshgrid(ds['lon'].values, ds['lat'].values)
    ds_seasons = ds.groupby('time.season')

    mu_season = {}
    sigma_season = {}
    ampl_season = {}
    phase_season = {}
    ampl_cov_season = {}
    phase_cov_season = {}
    average_cycle_season = {}

    for season_i, season_ds in ds_seasons:
#         if season_i != 'DJF':
#             break
        print(season_i)
        lst_da = compute_lst_array(season_ds, 
                                   bin_interval = grid_time_resolution_hours, # time_resolution_hours,
                                   bin_bool = True, 
                                   round_precision = round_precision,
                                   lon_mesh = lon_mesh, 
                                   lat_mesh = lat_mesh,
                                   field_id = field_id)
        lst_array = lst_da.values
        lst_array = lst_array.astype(np.float32)
        del lst_da
        lst_array[lst_array == 0] = 24
        
        field_season_array = season_ds[field_id].values
        field_season_array =  field_season_array.astype(np.float32)
#         field_season_mean = field_season_array.mean(axis = 0)
        print('float32')
        # compute mu_ij
#         mu_ij = np.zeros(field_season_array.shape[-2:])
        mu_ij = np.zeros(field_season_array.shape[-2:])
#         sigma_ij = np.zeros((len(hour_bins),) + field_season_array.shape[-2:])
        
        f_bar_ks = {}
        for hour_i in grid_hour_bins:
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_k[np.isnan(f_bar_k)] = 0
            f_bar_ks[hour_i] = f_bar_k
            mu_ij  += f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))
        average_cycle_season[season_i] = hour_means


        mu_ij = (1/len(grid_hour_bins))*mu_ij
#         mu_ij = (1/8)*mu_ij
        mu_season[season_i] = mu_ij
    
        sigma_ij = np.zeros(field_season_array.shape[-2:])
        for hour_i in grid_hour_bins:
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_k[np.isnan(f_bar_k)] = 0
            sigma_ij += np.square(f_bar_k - mu_ij)

        sigma_ij = np.sqrt((1/(len(grid_hour_bins) - 1))*sigma_ij)
#         sigma_ij = np.sqrt((1/7)*sigma_ij)
#         sigma_ij = np.nanstd(sigma_ij, axis = 0)

        sigma_season[season_i] = sigma_ij
    
        lst_da = compute_lst_array(season_ds, 
                                   bin_interval = time_resolution_hours,
                                   bin_bool = True, 
                                   round_precision = round_precision,
                                   lon_mesh = lon_mesh, 
                                   lat_mesh = lat_mesh,
                                   field_id = field_id)
        lst_array = lst_da.values
        del lst_da
        lst_array[lst_array == 0] = 24
        field_season_array = season_ds[field_id].values

#         field_season_mean = field_season_array.mean(axis = 0)

        # compute hourly grid means needed for cos fit
        f_bar_ks = {}
        for ii in range(len(hour_bins)):
            hour_i = hour_bins[ii]
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_ks[hour_i] = f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))

        print('Performing Cos Fit')
#         res = cos_fit_grid_average(hour_means, hour_bins)
#         pickle_path = '/export/data1/cchristo/diurnal_analysis_results/CMIP6_bin_pt2_hr_round_closest_gpm_overlap_2_modes_test/scrap/mm.pickle'
#         with open(pickle_path, 'wb') as f:
#             pickle.dump([hour_means, hour_bins], f)
            
        res = cos_fit_2_mode_grid_average(hour_means, hour_bins)
        print('Finished Cos Fit')
        ampl_season[season_i], phase_season[season_i] = res[0], res[1]
        ampl_cov_season[season_i], phase_cov_season[season_i] = res[2], res[3]
        
    out_ds = xr.Dataset()
    out_ds['mu_season'] = make_da_from_dict(mu_season, ds)
    out_ds['sigma_season'] = make_da_from_dict(sigma_season,ds)
    out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds)
    out_ds['phase_season'] = make_da_from_dict(phase_season,ds)
    out_ds['ampl_cov_season'] = make_da_from_dict(ampl_cov_season, ds)
    out_ds['phase_cov_season'] = make_da_from_dict(phase_cov_season,ds)

    out_ds_means = xr.Dataset()
    out_ds_means = make_da_from_dict_time(average_cycle_season, ds, grid_hour_bins)
    out_ds_means = out_ds_means.to_dataset(name = field_id + '_mean')
       

    return (out_ds, out_ds_means)
    
##### Do full analysis 

def diurnal_analysis(ds, 
                     field_id, 
                     grid_time_resolution_hours = 3, 
                     time_resolution_hours = 1,
                     round_precision = 2):
    ''' Assumes ordering (time, lat, lon), with those field/coordinate names.
    Given dataset compute mean, standard deviation, amplitude, and phase. '''
    hour_bins = np.arange(time_resolution_hours, 24 + time_resolution_hours, time_resolution_hours)
    hour_bins = np.round(hour_bins, round_precision)
    grid_hour_bins = np.arange(grid_time_resolution_hours, 24 + grid_time_resolution_hours, grid_time_resolution_hours)
    grid_hour_bins = np.round(grid_hour_bins, round_precision)
    
    lon_mesh, lat_mesh = np.meshgrid(ds['lon'].values, ds['lat'].values)
    ds_seasons = ds.groupby('time.season')

    mu_season = {}
    sigma_season = {}
    ampl_season = {}
    phase_season = {}
    ampl_cov_season = {}
    phase_cov_season = {}
    average_cycle_season = {}

    for season_i, season_ds in ds_seasons:
#         if season_i != 'DJF':
#             break
        print(season_i)
        lst_da = compute_lst_array(season_ds, 
                                   bin_interval = grid_time_resolution_hours, # time_resolution_hours,
                                   bin_bool = True, 
                                   round_precision = round_precision,
                                   lon_mesh = lon_mesh, 
                                   lat_mesh = lat_mesh,
                                   field_id = field_id)
        lst_array = lst_da.values
        lst_array = lst_array.astype(np.float32)
        del lst_da
        lst_array[lst_array == 0] = 24
        
        field_season_array = season_ds[field_id].values
        field_season_array =  field_season_array.astype(np.float32)
#         field_season_mean = field_season_array.mean(axis = 0)
        print('float32')
        # compute mu_ij
#         mu_ij = np.zeros(field_season_array.shape[-2:])
        mu_ij = np.zeros(field_season_array.shape[-2:])
#         sigma_ij = np.zeros((len(hour_bins),) + field_season_array.shape[-2:])
        
        f_bar_ks = {}
        for hour_i in grid_hour_bins:
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_k[np.isnan(f_bar_k)] = 0
            f_bar_ks[hour_i] = f_bar_k
            mu_ij  += f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))
        average_cycle_season[season_i] = hour_means


        mu_ij = (1/len(grid_hour_bins))*mu_ij
#         mu_ij = (1/8)*mu_ij
        mu_season[season_i] = mu_ij
    
        sigma_ij = np.zeros(field_season_array.shape[-2:])
        for hour_i in grid_hour_bins:
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_k[np.isnan(f_bar_k)] = 0
            sigma_ij += np.square(f_bar_k - mu_ij)

        sigma_ij = np.sqrt((1/(len(grid_hour_bins) - 1))*sigma_ij)
#         sigma_ij = np.sqrt((1/7)*sigma_ij)
#         sigma_ij = np.nanstd(sigma_ij, axis = 0)

        sigma_season[season_i] = sigma_ij
    
        lst_da = compute_lst_array(season_ds, 
                                   bin_interval = time_resolution_hours,
                                   bin_bool = True, 
                                   round_precision = round_precision,
                                   lon_mesh = lon_mesh, 
                                   lat_mesh = lat_mesh,
                                   field_id = field_id)
        lst_array = lst_da.values
        del lst_da
        lst_array[lst_array == 0] = 24
        field_season_array = season_ds[field_id].values

#         field_season_mean = field_season_array.mean(axis = 0)

        # compute hourly grid means needed for cos fit
        f_bar_ks = {}
        for ii in range(len(hour_bins)):
            hour_i = hour_bins[ii]
            masked_field = np.where(lst_array == hour_i, field_season_array, np.nan)

            # mean for a given season, LST
            f_bar_k = np.nanmean(masked_field, axis = 0)
            f_bar_ks[hour_i] = f_bar_k

        hour_means = np.stack(list(f_bar_ks.values()))

        print('Performing Cos Fit')
#         res = cos_fit_grid(field_season_array, lst_array, hour_means, hour_bins)
        res = cos_fit_grid_average(hour_means, hour_bins)
        print('Finished Cos Fit')
        ampl_season[season_i], phase_season[season_i] = res[0], res[1]
        ampl_cov_season[season_i], phase_cov_season[season_i] = res[2], res[3]
        
    out_ds = xr.Dataset()
    out_ds['mu_season'] = make_da_from_dict(mu_season, ds)
    out_ds['sigma_season'] = make_da_from_dict(sigma_season,ds)
    out_ds['ampl_season'] = make_da_from_dict(ampl_season, ds)
    out_ds['phase_season'] = make_da_from_dict(phase_season,ds)
    out_ds['ampl_cov_season'] = make_da_from_dict(ampl_cov_season, ds)
    out_ds['phase_cov_season'] = make_da_from_dict(phase_cov_season,ds)

    out_ds_means = xr.Dataset()
    out_ds_means = make_da_from_dict_time(average_cycle_season, ds, grid_hour_bins)
    out_ds_means = out_ds_means.to_dataset(name = field_id + '_mean')
       

    return (out_ds, out_ds_means)
    
    
    
    
##### Utility functions for plotting results 

def make_four_panel(season_ds, 
                    title = r'$\Phi$',
                    axis = None,
                    cmap = plt.get_cmap('twilight'),
                    vmin = None,
                    vmax = None,
                    save_fig_path = None):
    ''' Make 4 panel plot. '''
    lats = season_ds['lat'].values
    lons = season_ds['lon'].values
    
    def set_axis(ax, axis_list = axis):
        ax.set_xlim(axis_list[:2])
        ax.set_ylim(axis_list[2:])
        

    fig, axes = plt.subplots(nrows=4, ncols=1, sharex = True, figsize = (8,10))
    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    
    for ii in range(len(axes.flat)):
        if seasons[ii] in season_ds['season']:
            m = Basemap(ax=axes.flat[ii], llcrnrlon = lons[0] , 
                llcrnrlat = lats[0], urcrnrlon = lons[-1], urcrnrlat = lats[-1])
            im = m.pcolormesh(lons, lats, season_ds.sel(season = seasons[ii]).values, cmap = cmap, vmin = vmin, vmax = vmax)

            m.drawcoastlines()

            m.drawparallels(range(int(np.floor(np.min(lats))), int(np.ceil(np.max(lats))), 30), 
                            labels=[0,1,1,1])
            m.drawmeridians(range(int(np.min(lons)), int(np.max(lons)), 40), 
                            labels=[0,0,0,1])
            # plt.colorbar()
            if ii == 0:
                axes.flat[ii].set_title(title)
            axes.flat[ii].set_ylabel(seasons[ii])
            if axis: 
                set_axis(axes.flat[ii])
        

    fig.colorbar(im, ax=axes.ravel().tolist())
#     plt.tight_layout()
    if save_fig_path:
        plt.savefig(save_fig_path, dpi = 300)
        
        
def make_single_plot(data,
                lats = None,
                lons = None, 
                ax = None,
                fig = None,
                title = r'$\Phi$',
                title_fontsize = 13,
                ylabel = None,
                subplot_label = None,
                xticks_bool = True,
                axis = None,
                cmap = plt.get_cmap('twilight'),
                cbar = None,
                cbar_frac = 0.01,
                cbar_ylabel = None,
                vmin = None,
                figsize = (12,4),
                vmax = None,
                save_fig_path = None):
    ''' Make single geographic plot. '''
    
    if isinstance(data, xr.DataArray):
        lats = data['lat'].values
        lons = data['lon'].values
        data = data.values
    
    if (fig is None) & (ax is None):
        fig = plt.figure(figsize = figsize)
        ax = plt.gca()
#     fig = plt.figure(figsize = (10,3))
    
    def set_axis(ax, axis_list = axis):
        ax.set_xlim(axis_list[:2])
        ax.set_ylim(axis_list[2:])
        

    m = Basemap(ax=ax, llcrnrlon = lons[0] , 
        llcrnrlat = lats[0], urcrnrlon = lons[-1], urcrnrlat = lats[-1])
    im = m.pcolormesh(lons, lats, data, cmap = cmap, vmin = vmin, vmax = vmax)

    m.drawcoastlines()
#     int(np.floor(np.min(lats)))
    if xticks_bool:
        parallel_labels = [0,0,0,1]
    else:
        parallel_labels = [0,0,0,0]
    m.drawparallels(range(-30, int(np.ceil(np.max(lats))), 30), 
                    labels=[0,1,0,1])
    m.drawmeridians(range(int(np.min(lons)), int(np.max(lons)), 30), 
                    labels=parallel_labels, rotation = 45)
    # plt.colorbar()

#     ax.set_ylabel('Latitude')
#     ax.set_xlabel('Longitude')
#     plt.xlabel('Longitude')
    plt.xticks(rotation = 45)
    if title:
        ax.set_title(title, 
                     weight = 'bold',
                     fontsize = title_fontsize)
    if not (ylabel is None):
        ax.set_ylabel(ylabel, weight = 'bold')
    if subplot_label: 
            ax.text(-8, 53, subplot_label, fontweight = 'bold')
    if axis: 
        set_axis(ax)
        
    if cbar:
        if type(cbar) == bool:
            cbar_obj = fig.colorbar(im, 
                     ax=ax, 
                     fraction = cbar_frac, 
                     pad = 0.088)    #pad = 0.09
        else:
            cbar_obj = fig.colorbar(cbar, 
                         ax=ax, 
                         fraction = cbar_frac, 
                         pad = 0.088)    #pad = 0.09
        cbar_obj.ax.set_xlabel(cbar_ylabel, fontweight = 'bold')
        
#     if cbar:
#         fig.colorbar(cbar, ax=ax, fraction = 0.01, pad = 0.09)
#     else:
#         fig.colorbar(im, ax=ax, fraction = 0.01, pad = 0.09)
        
    plt.tight_layout()

    
    if save_fig_path:
        plt.savefig(save_fig_path, dpi = 300)
    
    
def make_da_from_dict(time_dict, ds, lat_id = 'lat', lon_id = 'lon', time_id = 'season'):
    ''' Given dict with arrays (no time dim) for each season, create a dataarray with the data. '''
    return xr.DataArray(np.stack(list(time_dict.values())), dims = (time_id, lat_id, lon_id), 
                        coords = {time_id:list(time_dict.keys()),
                        lon_id: ds[lon_id].values, 
                        lat_id: ds[lat_id].values})

def make_da_from_dict_full(time_dict, ds,lat_id = 'lat', lon_id = 'lon', time_id = 'season'):
    ''' Given dict with arrays (no time dim) for each season, create a dataarray with the data. '''
    return xr.DataArray(np.stack(list(time_dict.values())), dims = (time_id, lat_id, lon_id), 
                        coords = {time_id:list(time_dict.keys()),
                        lon_id: ds[lon_id].values, 
                        lat_id: ds[lat_id].values})

def make_da_from_dict_time(time_dict, ds, hour_bins, time_id = 'time', lat_id = 'lat', lon_id = 'lon', season_id = 'season'):
    ''' Given dict with arrays for each season, create a dataarray with the data. '''
    return xr.DataArray(np.stack(list(time_dict.values())), dims = (season_id,
                                                                      time_id,
                                                                      lat_id, 
                                                                      lon_id), 
                        coords = {season_id:list(time_dict.keys()),
                        time_id: hour_bins,
                        lon_id: ds[lon_id].values, 
                        lat_id: ds[lat_id].values})



def stationary_bootstrap_int(datasize, samplesize, blocksize, samples_int = 48, nboot = 200):
    
    """
    Create stationary bootstrap with integer number of days. 
    bootstrap_indices = stationary_bootstrap(datasize, samplesize, blocksize, nboot)
    samples_int - samples corresponding to int. (ie if 24 samples in 1 day, samples_int = 24)
    As in Poltis and Romano 1984: Returns `nboot` x `samplesize` array of indices from 
    0 to `datasize` using blocks of size given by geometric distribution of size `blocksize`.
    """
    bootstrap_indices = np.zeros((nboot, samplesize), dtype=int)
    for i in range(nboot):
        # size of sample follows geometric distribution with mean nsamples_mean
        blocksizes = np.random.geometric(1/blocksize, size=datasize) 
        blocksizes[blocksizes < samples_int] = samples_int
        blocksizes = np.round(blocksizes/samples_int).astype(int)
        blocksizes = (blocksizes * samples_int).astype(int)
#         return blocksizes
        # random starting indices of each block
        indices = np.random.randint(datasize, size=datasize)
        # construct pseudodata
        j = 0
        k = 0
        while j < samplesize:
            I = indices[k]
            L = blocksizes[k]
            for n in range(L):
                if j + n == samplesize:
                    break
                bootstrap_indices[i, j + n] = (I + n) % datasize
            j += L
            k += 1
    return bootstrap_indices


def stationary_bootstrap(datasize, samplesize, blocksize, nboot = 200):
    
    """
    bootstrap_indices = stationary_bootstrap(datasize, samplesize, blocksize, nboot)
    As in Poltis and Romano 1984: Returns `nboot` x `samplesize` array of indices from 
    0 to `datasize` using blocks of size given by geometric distribution of size `blocksize`.
    """
    bootstrap_indices = np.zeros((nboot, samplesize), dtype=int)
    for i in range(nboot):
        # size of sample follows geometric distribution with mean nsamples_mean
        blocksizes = np.random.geometric(1/blocksize, size=datasize) 
        # random starting indices of each block
        indices = np.random.randint(datasize, size=datasize)
        # construct pseudodata
        j = 0
        k = 0
        while j < samplesize:
            I = indices[k]
            L = blocksizes[k]
            for n in range(L):
                if j + n == samplesize:
                    break
                bootstrap_indices[i, j + n] = (I + n) % datasize
            j += L
            k += 1
    return bootstrap_indices
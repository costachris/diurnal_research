'''Utilities for loading specific models downloaded with CMIPX_downloader. 
Assumes all desired models/ensemble members are in the same directory'''

import os
# from tqdm import tqdm



def get_unique_models(path_to_model_dirs):
    '''Given path to dir containing model files, get unique names of models. '''
    model_list = os.listdir(path_to_model_dirs)
    desired_model_fnames = []
    for model_i in range(len(model_list)):
        filename_split = model_list[model_i].split('_')
        desired_model_fnames.append(filename_split[2])
    return set(desired_model_fnames)


# desired_model_fnames


def get_path_to_desired_model_cmip6(path_to_model_dirs,
                              desired_model = 'GFDL-CM4',
                              desired_ensemble_member = 'r1i1p1f1',
                              desired_grid_types = ('gn', 'gr', 'gr1')): 
    '''Get list of model file paths matching criteria for cmip6.
    Args:
        path_to_model_dirs - str
            path to dir containing all models, including trailing backslash
            
        desired_grid_types - tuple
            desired grid types in decreasing order of preference. CMIP6 has different grid types, 
            CMIP5 has 1 type.
    '''
    model_list = os.listdir(path_to_model_dirs)
    desired_model_fnames = []
    for model_i in range(len(model_list)):
        filename_split = model_list[model_i].split('_')

        if (filename_split[-3] == desired_ensemble_member) & \
           (filename_split[2] == desired_model):
            for grid_type_i in desired_grid_types:
                if filename_split[-2] == grid_type_i:
                    desired_model_fnames.append(path_to_model_dirs + model_list[model_i])
                else:
                    continue
    return sorted(desired_model_fnames)

def get_path_to_desired_model_cmip5(path_to_model_dirs,
                              desired_model = 'GFDL-CM4',
                              desired_ensemble_member = 'r1i1p1',): 
    '''Get list of model file paths matching criteria for cmip5.
    Args:
        path_to_model_dirs - str
            path to dir containing all models, including trailing backslash
            
        desired_grid_types - tuple
            desired grid types in decreasing order of preference. CMIP6 has different grid types, 
            CMIP5 has 1 type.
    '''
    model_list = os.listdir(path_to_model_dirs)
    desired_model_fnames = []
    for model_i in range(len(model_list)):
        filename_split = model_list[model_i].split('_')

        if (filename_split[-2] == desired_ensemble_member) & \
           (filename_split[2] == desired_model):
            desired_model_fnames.append(path_to_model_dirs + model_list[model_i])


    return sorted(desired_model_fnames)
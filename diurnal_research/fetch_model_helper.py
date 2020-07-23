'''Utilities for loading specific models downloaded with CMIPX_downloader. 
Assumes all desired models/ensemble members are in the same directory'''

import os
# from tqdm import tqdm



def get_unique_models(path_to_model_dirs):
    '''Given path to dir containing model files, get unique names of models. 
    Returns
    --------
        desired_model_fnames - list
            unique CMIP models in dir
    
    '''
    model_list = os.listdir(path_to_model_dirs)
    desired_model_fnames = []
    for model_i in range(len(model_list)):
        filename_split = model_list[model_i].split('_')
        desired_model_fnames.append(filename_split[2])
    return set(desired_model_fnames)


def _get_ensemble_from_path_cmip6(path):
    return path.split('_')[-3]

def _get_ensemble_from_path_cmip5(path):
    return path.split('_')[-2]

# desired_model_fnames
def get_unique_ensembles_cmip6(path_to_model_dirs):
    ''' Given path to CMIP6 models, get unique ensembles for each model 
    
    Returns
    --------
        unqiue_ensembles - dict
            dictionary mapping unique models to available ensemble members
    '''
    unique_model_names = get_unique_models(path_to_model_dirs)
    model_to_ensemble_dict = {}
    for model_ii in unique_model_names:
        paths_for_model = get_path_to_desired_model_cmip6(path_to_model_dirs,
                                                          model_ii,
                                                          desired_ensemble_member = None)
        ensembles = [_get_ensemble_from_path_cmip6(path) for path in paths_for_model]
        model_to_ensemble_dict[model_ii] = set(ensembles)
        
    return model_to_ensemble_dict


def get_path_to_desired_model_cmip6(path_to_model_dirs,
                              desired_model = 'GFDL-CM4',
                              desired_ensemble_member = ('r1i1p1f1', 'r1i1p2f1','r1i1p1f2'),
                              desired_grid_types = ('gn', 'gr', 'gr1', 'gr2')): 
    '''Get list of model file paths matching criteria for cmip6.
    Args:
        path_to_model_dirs - str
            path to dir containing all models, including trailing backslash
            
        desired_ensemble_member - tuple 
                desired ensemble members in decreasing order of preference
            
        desired_grid_types - tuple
            desired grid types in decreasing order of preference. CMIP6 has different grid types, 
            CMIP5 has 1 type.
    '''
    model_list = os.listdir(path_to_model_dirs)
    desired_model_fnames = []
    for model_i in range(len(model_list)):
        filename_split = model_list[model_i].split('_')
        # if an ensemble member is requested, return only that ensemble member
        if desired_ensemble_member:
            for desired_ensemble_member_ii in desired_ensemble_member:
                if (filename_split[-3] == desired_ensemble_member_ii) & \
                   (filename_split[2] == desired_model):
                    for grid_type_i in desired_grid_types:
                        if filename_split[-2] == grid_type_i:
                            desired_model_fnames.append(path_to_model_dirs + model_list[model_i])
                            break
                        else:
                            continue
                pass
        
        # if an ensemble member not requested, return all ensemble members
        else: 
            if (filename_split[2] == desired_model):
                for grid_type_i in desired_grid_types:
                    if filename_split[-2] == grid_type_i:
                        desired_model_fnames.append(path_to_model_dirs + model_list[model_i])
                        break
                    else:
                        continue
    # if dir contains multiple grid types, we only want to select the most desired (otherwise we'll get repeated files)
    desired_model_fnames_final = []
    grid_types_in_desired_paths = [f_path.split('_')[-2] for f_path in desired_model_fnames]
    for grid_type_ii in desired_grid_types:
        if grid_type_ii in grid_types_in_desired_paths:
            for i in range(len(grid_types_in_desired_paths)):
                if grid_types_in_desired_paths[i] == grid_type_ii:
                    desired_model_fnames_final.append(desired_model_fnames[i])
            break
            
    # if dir contains multiple esemble members, we only want to select the most desired (otherwise we'll get repeated files)
    all_desired_file_paths = []
    ensembles_in_desired_paths = [f_path.split('_')[-3] for f_path in desired_model_fnames_final]
    for ensemble_member_ii in desired_ensemble_member:
        if ensemble_member_ii in ensembles_in_desired_paths:
            for i in range(len(ensembles_in_desired_paths)):
                if ensembles_in_desired_paths[i] == ensemble_member_ii:
                    all_desired_file_paths.append(desired_model_fnames_final[i])
            break
            
    return sorted(all_desired_file_paths)


def get_unique_ensembles_cmip5(path_to_model_dirs):
    ''' Given path to CMIP5 models, get unique ensembles for each model 
    
    Returns
    --------
        unqiue_ensembles - dict
            dictionary mapping unique models to available ensemble members
    '''
    unique_model_names = get_unique_models(path_to_model_dirs)
    model_to_ensemble_dict = {}
    for model_ii in unique_model_names:
        paths_for_model = get_path_to_desired_model_cmip5(path_to_model_dirs,
                                                          model_ii,
                                                          desired_ensemble_member = None)
        ensembles = [_get_ensemble_from_path_cmip5(path) for path in paths_for_model]
        model_to_ensemble_dict[model_ii] = set(ensembles)
        
    return model_to_ensemble_dict




def get_path_to_desired_model_cmip5(path_to_model_dirs,
                              desired_model = 'GFDL-CM3',
                              desired_ensemble_member = ('r1i1p1','r6i1p1','r2i1p1')): 
    '''Get list of model file paths matching criteria for cmip5.
    Args:
        path_to_model_dirs - str
            path to dir containing all models, including trailing backslash
        
        desired_ensemble_member - tuple 
                desired ensemble members in decreasing order of preference
            
        desired_grid_types - tuple
            desired grid types in decreasing order of preference. CMIP6 has different grid types, 
            CMIP5 has 1 type.
    '''
    model_list = os.listdir(path_to_model_dirs)
    desired_model_fnames = []
    for model_i in range(len(model_list)):
        filename_split = model_list[model_i].split('_')
        # if an ensemble member is requested, return only that ensemble member
        if desired_ensemble_member:
            for desired_ensemble_member_ii in desired_ensemble_member:
                if (filename_split[-2] == desired_ensemble_member_ii) & \
                   (filename_split[2] == desired_model):
                    desired_model_fnames.append(path_to_model_dirs + model_list[model_i])
                pass
         # if an ensemble member not requested, return all ensemble members
        else: 
            if (filename_split[2] == desired_model):
                desired_model_fnames.append(path_to_model_dirs + model_list[model_i])

    return sorted(desired_model_fnames)
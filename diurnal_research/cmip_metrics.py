'''Externally computed metrics for CMIP models, by name (ie. ECS, TCS, error) 
'Context for interpreting equilibrium climate sensitivity and transient climate response from the CMIP6 Earth system models' Science, 2020
https://advances.sciencemag.org/content/6/26/eaba1981 '''
import pandas as pd
###### CMIP6 metrics


cmip6_ecs = {
 'ACCESS-CM2': 4.7,
 'ACCESS-ESM1-5': 3.9, 
 'AWI-CM-1-1-MR': 3.2,
#  'AWI-ESM-1-1-LR':  # model not in paper
 'BCC-CSM2-MR': 3.0,
 'CanESM5': 5.6,
 #  'CMCC-CM2-SR5':  # model not in paper
 'CNRM-CM6-1': 4.8,
 'CNRM-CM6-1-HR': 4.3,
 'CNRM-ESM2-1': 4.8,
 'EC-Earth3': 4.3,
 'EC-Earth3-Veg': 4.3,
#  'EC-Earth3-Veg-LR': 4.3, # model not in paper
#  'FGOALS-g3': 3.0,  # name doesn't quite match
 
 'GFDL-CM4': 3.9, 
 'GFDL-ESM4': 2.6, 
 'GISS-E2-1-G': 2.7, 
 'IPSL-CM6A-LR': 4.6,
 'KACE-1-0-G': 4.5,
 'MIROC6': 2.6,
 'MIROC-ES2L': 2.7,
#   'MPI-ESM-1-2-HAM' # model not in paper
 'MPI-ESM1-2-HR': 3.0,
 'MPI-ESM1-2-LR': 3.0, 
 'MRI-ESM2-0': 3.2,
 'NESM3': 4.7,
 'SAM0-UNICON': 3.7,
 'UKESM1-0-LL': 5.3,
}


cmip6_tcr = {
 'ACCESS-CM2': 2.1,
 'ACCESS-ESM1-5': 2.0, 
 'AWI-CM-1-1-MR': 2.0,
 # 'AWI-ESM-1-1-LR': # model not in paper
 'BCC-CSM2-MR': 1.7,
 'CanESM5': 2.7,
 #  'CMCC-CM2-SR5': 
 'CNRM-CM6-1': 2.1,
 'CNRM-CM6-1-HR': 2.5,
 'CNRM-ESM2-1': 1.9,
 # 'EC-Earth3': , # tcr not reported in paper
 'EC-Earth3-Veg': 2.6,
 # 'EC-Earth3-Veg-LR': , # model not in paper
 # 'FGOALS-g3': 2.1, # name doesn't quite match
 'GFDL-CM4': 2.1, 
 'GFDL-ESM4': 1.6, 
 'GISS-E2-1-G': 1.8, 
 'IPSL-CM6A-LR': 2.3,
 'KACE-1-0-G': 1.4,
 'MIROC6': 1.6,
 'MIROC-ES2L': 1.6,
#   'MPI-ESM-1-2-HAM'  # model not in paper
 'MPI-ESM1-2-HR': 1.7,
 'MPI-ESM1-2-LR': 1.8, 
 'MRI-ESM2-0': 1.6,
 'NESM3': 2.7,
 'SAM0-UNICON': 2.3,
 'UKESM1-0-LL': 2.8,
}



###### CMIP5 metrics


cmip5_ecs = {
 'ACCESS1-0': 3.8,
 'ACCESS1-3': 3.5,
#  'CMCC-CM': ,
 'CNRM-CM5': 3.3,
 'FGOALS-g2': 3.4,
 'FGOALS-s2': 4.2,
 'GFDL-CM3': 4.0,
 'GFDL-ESM2G': 2.4,
 'GFDL-ESM2M': 2.4,
 'GISS-E2-H': 2.3,
 'GISS-E2-R': 2.1,
 'HadGEM2-ES': 4.6,
 'IPSL-CM5A-LR': 4.1,
#  'IPSL-CM5A-MR':,
 'MIROC-ESM': 4.7,
#  'MIROC-ESM-CHEM',
#  'MIROC4h',
 'MIROC5': 2.7,
 'MRI-CGCM3': 2.6,
#  'MRI-ESM1': ,
 'NorESM1-M': 2.8,
 'inmcm4': 2.1,
}


cmip5_tcr = {
 'ACCESS1-0': 1.9,
 'ACCESS1-3': 1.6,
#  'CMCC-CM': ,
 'CNRM-CM5': 2.0,
 'FGOALS-g2': 1.4,
 'FGOALS-s2': 2.4,
 'GFDL-CM3': 1.9,
 'GFDL-ESM2G': 1.1,
 'GFDL-ESM2M': 1.4,
 'GISS-E2-H': 1.7,
 'GISS-E2-R': 1.5,
 'HadGEM2-ES': 2.5,
 'IPSL-CM5A-LR': 2.0,
 'IPSL-CM5A-MR': 2.0,
 'MIROC-ESM': 2.2,
#  'MIROC-ESM-CHEM',
#  'MIROC4h',
 'MIROC5': 1.4,
 'MRI-CGCM3': 1.6,
#  'MRI-ESM1': , # name doesn't quite match table
 'NorESM1-M': 1.4,
 'inmcm4': 1.3,
}



cmip6_ecs_df = pd.DataFrame.from_dict(cmip6_ecs, orient = 'index')
cmip6_ecs_df.columns = ['ECS',]

cmip6_tcr_df = pd.DataFrame.from_dict(cmip6_tcr, orient = 'index')
cmip6_tcr_df.columns = ['TCR',]

cmip6_sensitivities = pd.merge(cmip6_ecs_df, cmip6_tcr_df, how = 'outer', left_index = True, right_index = True)

cmip5_ecs_df = pd.DataFrame.from_dict(cmip5_ecs, orient = 'index')
cmip5_ecs_df.columns = ['ECS',]

cmip5_tcr_df = pd.DataFrame.from_dict(cmip5_tcr, orient = 'index')
cmip5_tcr_df.columns = ['TCR',]

cmip5_sensitivities = pd.merge(cmip5_ecs_df, cmip5_tcr_df, how = 'outer', left_index = True, right_index = True)
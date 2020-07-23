'''Configurations for diurnal analysis accross all CMIP 5/6 models. '''


# unique cmip models with 3hr pr data
unique_pr_cmip6_models = \
    ('ACCESS-CM2',
     'ACCESS-ESM1-5',
     'AWI-CM-1-1-MR',
     'AWI-ESM-1-1-LR',
     'BCC-CSM2-MR',
     'CMCC-CM2-SR5',
     'CNRM-CM6-1',
     'CNRM-CM6-1-HR',
     'CNRM-ESM2-1',
     'CanESM5',
     'EC-Earth3',
     'EC-Earth3-Veg',
     'EC-Earth3-Veg-LR',
     'FGOALS-g3',
     'GFDL-CM4',
     'GFDL-ESM4',
     'GISS-E2-1-G',
     'IPSL-CM6A-LR',
     'KACE-1-0-G',
     'MIROC-ES2L',
     'MIROC6',
     'MPI-ESM-1-2-HAM',
     'MPI-ESM1-2-HR',
     'MPI-ESM1-2-LR',
     'MRI-ESM2-0',
     'NESM3',
     'SAM0-UNICON',
     'UKESM1-0-LL')



cmip6_to_cmip5_map = {
    'GFDL-ESM4':'GFDL-ESM2M',
    'GFDL-CM4':'GFDL-CM3',
    'FGOALS-g3':'FGOALS-g2',
    'IPSL-CM6A-LR':'IPSL-CM5A-LR',
    'MIROC6':'MIROC5',
    'MRI-ESM2-0': 'MRI-ESM1',
    
}
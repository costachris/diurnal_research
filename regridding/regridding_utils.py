import os
from glob import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import xesmf as xe
import logging


regrid_weights_dir = '/export/data1/cchristo/regridding_weights/'


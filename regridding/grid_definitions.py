'''Grid definitions used for regridding. '''

import numpy as np
import xarray as xr

grids = {
    'grid1': xr.Dataset({'lat': (['lat'], np.arange(-60, 60, 0.5)),
                     'lon': (['lon'], np.arange(0, 360, 0.5)),})
}


'''Grid definitions used for regridding. '''

import numpy as np
import xarray as xr

grids = {
    'grid1': xr.Dataset({'lat': (['lat'], np.arange(-60, 60, 0.5)),
                     'lon': (['lon'], np.arange(0, 360, 0.5)),}),
    'ryan_1deg' :  xr.Dataset({'lat': (['lat'], np.arange(-89.5, 89.5 + 1, 1.0)),
                     'lon': (['lon'], np.arange(-179.5,179.5 + 1, 1.0)),})
}


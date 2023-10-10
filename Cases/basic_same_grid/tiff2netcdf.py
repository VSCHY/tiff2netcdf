import rioxarray
import numpy as np
import numpy.ma as ma
import xarray as xr

import time
# LOAD ALL TIFF FILE
dfile = "/home/aschrapffer/DATA/FLOOD/Floodmap_EFAS/floodmap_GL/floodMapGL_rp100y/floodMapGL_rp100y.tif"

t0=time.time()
d = rioxarray.open_rasterio(dfile)
data = d.values[0,:,:]

lon = d["x"].values
lat = d["y"].values

data[data<0] = 0

ds = xr.Dataset({
    "Flood100Y": xr.DataArray(
        data   = data,
        dims   = ['latitude',"longitude"],
        coords = {'latitude':lat,"longitude":lon},
        attrs  = {'units': "flood_depth_m"},
        )},
        attrs = {})

ENC = {"zlib": True, "complevel": 1, "fletcher32": True}
ds.to_netcdf("Flood_100Y.nc", format= 'NETCDF4', engine='h5netcdf', encoding = {"Flood100Y":ENC})





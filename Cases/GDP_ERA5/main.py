import sys
sys.path.append("../../")
from src import Grid, GDP_mask, ERA5_mask
import multiprocessing
import numpy as np
import time
import xarray as xr

######################

dlandsea = "/home/aschrapffer/DATA/ERA5/LANDSEAMASK/landsea.nc"
dgdp = "/home/aschrapffer/DATA/ECONOMIC/GDP_HDI_POP/GDP_PPP_1990_2015_5arcmin_v2.nc"

d_out_file = "/home/aschrapffer/Documents/tiff2netcdf/Cases/GDP_ERA5/Output/GDP_2015_ERA5.nc"

######################

gr_gdp = Grid(dgdp, extent = None)
gr_era5 = Grid(dlandsea, extent = None)

mask_gdp = GDP_mask(dgdp)
gr_gdp.load_mask(mask_gdp)
gr_gdp.select()

mask_era5 = ERA5_mask(dlandsea)
gr_era5.load_mask(mask_era5)
gr_era5.select()

####################################################

# In this case we add the GDP to the closest point being land: is there some case the point is too far ? 

lons = np.array([c[0] for c in gr_era5.centers])
lats = np.array([c[1] for c in gr_era5.centers])

def get_point(index):
    lon0, lat0 = gr_gdp.centers[index]
    out = (lons-lon0)**2 + (lats-lat0)**2
    i = np.argmin(out)
    return i

####################################################

nbcore = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=nbcore-1)
#
t0 = time.time()
P_list = pool.map(get_point, np.arange(gr_gdp.n_elements))
t1 = time.time()
#
print(f"{t1-t0:.2f}", "s.")

############
# Export output
data = gr_gdp.get_data("GDP_PPP")
data = data[-1,:,:] # Last year # note in variable name

data_out = np.zeros((gr_era5.nj,gr_era5.ni))
for index, p in enumerate(P_list):
    j,i = gr_era5.ind_land[p]
    jg,ig = gr_gdp.ind_land[index]
    data_out[j,i] += data[jg,ig]

ds = xr.Dataset({
        "GDP_2015": xr.DataArray(
            data   = data_out,
            dims   = ['latitude',"longitude"],
            coords = {'latitude':gr_era5.xarlat,"longitude":gr_era5.xarlon},
            attrs  = {"long_name":f"GDP PPP 2015", 'units': '-'},
            )},
        attrs = {})
ds.to_netcdf(d_out_file)
ds.close()
del ds; del data_out

#gr.create_output(P_list, "cropfrac", d_out_file)
t2 = time.time()
print(f"{t2-t1:.2f}", "s.")

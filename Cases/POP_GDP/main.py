import sys
sys.path.append("../../")
from src import Grid, POP_mask, GDP_mask
import multiprocessing
import numpy as np
import time
import xarray as xr

######################

dgdp = "/home/aschrapffer/DATA/ECONOMIC/GDP_HDI_POP/GDP_PPP_1990_2015_5arcmin_v2.nc"
dpop = "/home/aschrapffer/DATA/ECONOMIC/GPW/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11_totpop_2pt5_min_nc/POPULATION_2020.nc"

d_out_file = "/home/aschrapffer/Documents/tiff2netcdf/Cases/POP_GDP/Output/POP_2020_GDP_grid.nc"

######################

gr_pop = Grid(dpop, extent = None)
gr_gdp = Grid(dgdp, extent = None)

mask_pop = POP_mask(dpop)
gr_pop.load_mask(mask_pop)
gr_pop.select()

mask_gdp = GDP_mask(dgdp)
gr_gdp.load_mask(mask_gdp)
gr_gdp.select()

####################################################

lons = np.array([c[0] for c in gr_gdp.centers])
lats = np.array([c[1] for c in gr_gdp.centers])
indexes = np.array([ii for ii in range(len(gr_gdp.centers))])

def get_point(index):
    #
    lon0, lat0 = gr_pop.centers[index]

    condition = (np.abs(lons-lon0)<0.3)*(np.abs(lats-lat0)<0.3)

    if len(condition) == 0:
        condition = (np.abs(lons-lon0)<1)*(np.abs(lats-lat0)<1)
        if len(condition) == 0:
            condition = (np.abs(lons-lon0)<4)*(np.abs(lats-lat0)<4)

    usecondition = False
    if len(condition) > 0:
        if True in condition:
            usecondition = True 

    if usecondition > 0:
        selected_lon = lons[condition]
        selected_lat = lats[condition]
        selected_indexes = indexes[condition]

        out = (selected_lon-lon0)**2 + (selected_lat-lat0)**2
        try:
            i = np.argmin(out)
        except:
            print(out)
            print(index)
            print(lon0,lat0)
            print(condition)
            aaa
        return selected_indexes[i]
    else:
        return None

# argmin of an empty sequence

####################################################

print("Start")
nbcore = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=nbcore)
#
t0 = time.time()
P_list = pool.map(get_point, np.arange(gr_pop.n_elements))
t1 = time.time()
#
print(f"{t1-t0:.2f}", "s.")


print(gr_pop.n_elements)
############
# Export output
data = gr_pop.get_data("POPULATION 2020")

data_out = np.zeros((gr_gdp.nj,gr_gdp.ni))
for index, p in enumerate(P_list):
    if p is not None:
        j,i = gr_gdp.ind_land[p]
        jg,ig = gr_pop.ind_land[index]
        data_out[j,i] += data[jg,ig]
    else:
        pass

ds = xr.Dataset({
        "GDP_2015": xr.DataArray(
            data   = data_out,
            dims   = ['latitude',"longitude"],
            coords = {'latitude':gr_gdp.xarlat,"longitude":gr_gdp.xarlon},
            attrs  = {"long_name":f"Fraction of cropland", 'units': '-'},
            )},
        attrs = {})
ds.to_netcdf(d_out_file)
ds.close()
del ds; del data_out

#gr.create_output(P_list, "cropfrac", d_out_file)
t2 = time.time()
print(f"{t2-t1:.2f}", "s.")

import sys
sys.path.append("../../")
from src import Shapefile
from src import Grid, GDP_mask
import multiprocessing
import numpy as np
import time

######################

d_flopros = "/home/aschrapffer/DATA/FLOOD/River_Coast_Protections/Scussolini_etal_Suppl_info/FLOPROS_shp_V1/FLOPROS_shp_V1.shp"
dgdp = "/home/aschrapffer/DATA/ECONOMIC/GDP_HDI_POP/GDP_PPP_1990_2015_5arcmin_v2.nc"
TEMP = "/home/aschrapffer/Documents/tiff2netcdf/TEMP/"
d_out_file = "/home/aschrapffer/Documents/tiff2netcdf/Cases/FLOPROS_2_GDP/Output/FLOPROS_GDPgrid.nc"


######################

sh = Shapefile(d_flopros)

######################

gr_gdp = Grid(dgdp, extent = None)
mask_gdp = GDP_mask(dgdp)
gr_gdp.load_mask(mask_gdp)
gr_gdp.select()

VARN_LIST = ['DL_Min_Co','DL_Max_Co', 'PL_Min_Co', 'PL_Max_Co', 
        'MerL_Riv', 'DL_Min_Riv','DL_Max_Riv', 
        'PL_Min_Riv', 'PL_Max_Riv', 'ModL_Riv']


####################################################

def get_protection(index):

    lon0 = gr_gdp.centers[index][0]
    lat0 = gr_gdp.centers[index][1]
    
    gdf = sh.get_poly_from_lonlat(lon0, lat0)
    if gdf is None:
        return [-99 for varn in VARN_LIST]
    else:
       L_out = [gdf[varn].iloc[0] for varn in VARN_LIST]
       return L_out

####################################################

nbcore = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=nbcore-1)
#
t0 = time.time()
P_list = pool.map(get_protection, np.arange(gr_gdp.n_elements))


t1 = time.time()
#
print(f"{t1-t0:.2f}", "s.")


# Set variable after variables

def reconstruct(index_plist):
    data_out = np.full((gr_gdp.nj,gr_gdp.ni),-999)
    for index, p in enumerate(P_list):
        j,i = gr_gdp.ind_land[index]
        data_out[j,i] = p[index_plist]    
    return data_out
    
EXPORTED_DATASET = {}

for index_plist, varn in enumerate(VARN_LIST):
    data_out = reconstruct(index_plist)
    EXPORTED_DATASET[varn] = xr.DataArray(
            data   = data_out,
            dims   = ['latitude',"longitude"],
            coords = {'latitude':gr_gdp.xarlat,"longitude":gr_gdp.xarlon},
            attrs  = {"long_name":varn, 'units': 'year of return period'},
            )

ds = xr.Dataset(EXPORTED_DATASET,
        attrs = {})
ds.to_netcdf(d_out_file)
ds.close()
del ds; del data_out

t2 = time.time()
print(f"{t2-t1:.2f}", "s.")

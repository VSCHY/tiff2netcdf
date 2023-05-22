import sys
sys.path.append("../../")
from src import Shapefile
from src import Raster
import multiprocessing
import numpy as np
import time

######################

d_aqueduct = "/home/aschrapffer/DATA/AQUEDUCT/Y2019M07D12_Aqueduct30_V01/baseline/annual/y2019m07d11_aqueduct30_annual_v01.gpkg"
dagri = "/home/aschrapffer/DATA/AGRICULTURE/Global_Cropland/Global_cropland_3km_2019.tif"
TEMP = "/home/aschrapffer/Documents/tiff2netcdf/TEMP/"
d_out_file = "/home/aschrapffer/Documents/tiff2netcdf/Cases/GlobalCropland_Aqueduct/Output/Global_Cropland_Aqueduct.gpkg"


######################

sh = Shapefile(d_aqueduct)

######################

rast = Raster(dagri, TEMP)

####################################################

def get_fraction(index):
    poly = sh.gdf.iloc[index]["geometry"]
    name = f"CropLand_{index}"
    out = rast.get_fraction(poly, name)
    return out

####################################################

nbcore = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=nbcore-1)
#
t0 = time.time()
P_list = pool.map(get_fraction, np.arange(sh.n_elements))
t1 = time.time()
#
print(f"{t1-t0:.2f}", "s.")
#
out_gdf = sh.gdf.copy()
out_gdf["Cropland fraction"] = P_list
out_gdf.to_file(d_out_file, driver="GPKG")


t2 = time.time()
print(f"{t2-t1:.2f}", "s.")

import sys
sys.path.append("../../")
from src import Grid, GDP_mask
from src import Raster
import multiprocessing
import numpy as np
import time

######################

dGDP = "/home/aschrapffer/DATA/ECONOMIC/GDP_HDI_POP/GDP_PPP_1990_2015_5arcmin_v2.nc"
dagri = "/home/aschrapffer/DATA/AGRICULTURE/Global_Cropland/Global_cropland_3km_2019.tif"
TEMP = "/home/aschrapffer/Documents/tiff2netcdf/TEMP/"
d_out_file = "/home/aschrapffer/Documents/tiff2netcdf/Cases/GlobalCropland_GDP/Output/Global_Cropland_GDP.nc"

######################

gr = Grid(dGDP, extent = None)

mask = GDP_mask(dGDP)
gr.load_mask(mask)
gr.select()

######################

rast = Raster(dagri, TEMP)

####################################################

def get_fraction(index):
    poly = gr.polylist[index]
    name = f"CropLand_{index}"
    out = rast.get_fraction(poly, name)
    return out 

####################################################

nbcore = multiprocessing.cpu_count()

pool = multiprocessing.Pool(processes=nbcore-1)
#
t0 = time.time()
P_list = pool.map(get_fraction, np.arange(gr.n_elements))
t1 = time.time()
#
print(f"{t1-t0:.2f}", "s.")
#
gr.create_output(P_list, "cropfrac", d_out_file)
t2 = time.time()
print(f"{t2-t1:.2f}", "s.")

import sys
sys.path.append("../../")
from src import Grid, ERA5_mask
from src import Raster
import multiprocessing
import numpy as np
import time

######################

dlandsea = "/home/aschrapffer/DATA/ERA5/LANDSEAMASK/landsea.nc"
dagri = "/home/aschrapffer/DATA/AGRICULTURE/Global_Cropland/Global_cropland_3km_2019.tif"
TEMP = "/home/aschrapffer/Documents/extrapolation_tiff/TEMP/"

######################

gr = Grid(dlandsea, extent = None)

mask = ERA5_mask(dlandsea)
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
gr.create_output(P_list, "cropfrac", "/home/aschrapffer/Documents/extrapolation_tiff/Cases/GlobalCropland/Output/Global_Cropland_ERA5.nc")
t2 = time.time()
print(f"{t2-t1:.2f}", "s.")

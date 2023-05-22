import sys
sys.path.append("../../")
from src import Grid, GDP_mask
from src import Raster
import multiprocessing
import numpy as np
import time
import subprocess

######################

dire = "/home/aschrapffer/Documents/tiff2netcdf/Cases/Aqueduct_Flood_100y_Riverine_Present_GDP/"

dGDP = "/home/aschrapffer/DATA/ECONOMIC/GDP_HDI_POP/GDP_PPP_1990_2015_5arcmin_v2.nc"
dflood = "/home/aschrapffer/DATA/FLOOD/AqueductFloods/Riverine/inunriver_historical_000000000WATCH_1980_rp00100.tif"
TEMP = "/home/aschrapffer/Documents/tiff2netcdf/TEMP/"
d_out_file = dire + "Output/Aqueduct_Flood_100y_Riverine_Present_GDP_gt0.nc"

######################

gr = Grid(dGDP, extent = None)

mask = GDP_mask(dGDP)
gr.load_mask(mask)
gr.select()

######################

rast = Raster(dflood, TEMP)

####################################################

def get_fraction(index):
    poly = gr.polylist[index]
    name = f"Coastal_100y_present_{index}"
    out = rast.get_fraction_threshold(poly, name, 0)
    return out 

####################################################

# RM LES LOG + recréer le dossier
block_size = 150000

n_part = gr.n_elements//block_size
nbcore = multiprocessing.cpu_count()

#
for i in range(n_part):
    print(i)
    pool = multiprocessing.Pool(processes=14)
    ind_land = np.arange(i*block_size, (i+1)*block_size)
    t0  = time.time()
    P_list = pool.map(get_fraction, ind_land)
    t1 = time.time()
    #
    print(f"{i+1}/{n_part+1} {t1-t0:.2f}", "s.")
    #
    gr.create_output(P_list, "flood_frac", d_out_file.replace(".nc",f"_{i}.nc"), ind_land = ind_land)
    pool.close()


print(n_part)
pool = multiprocessing.Pool(processes=14)
ind_land = np.arange(n_part*block_size, gr.n_elements)
t0  = time.time()
P_list = pool.map(get_fraction, ind_land)
t1 = time.time()
#
print(f"{i+2}/{n_part+1} {t1-t0:.2f}", "s.")
#
gr.create_output(P_list, "flood_frac", d_out_file.replace(".nc",f"_{n_part}.nc"), ind_land = ind_land)
pool.close()




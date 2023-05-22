import sys
sys.path.append("../../")
from src import Shapefile
from src import Grid
from src import GDP_mask
import multiprocessing
import numpy as np
import time
import subprocess

######################

d_aqueduct = "/home/aschrapffer/DATA/AQUEDUCT/Y2019M07D12_Aqueduct30_V01/baseline/annual/y2019m07d11_aqueduct30_annual_v01.gpkg"
dGDP = "/home/aschrapffer/DATA/ECONOMIC/GDP_HDI_POP/GDP_PPP_1990_2015_5arcmin_v2.nc"
TEMP = "/home/aschrapffer/Documents/tiff2netcdf/TEMP/"
d_out_file = "/home/aschrapffer/Documents/tiff2netcdf/Cases/GDP_Aqueduct/Output/GDP_Aqueduct.gpkg"
d_log = "/home/aschrapffer/Documents/tiff2netcdf/Cases/GDP_Aqueduct/Log/"

######################

sh = Shapefile(d_aqueduct)

######################

gr = Grid(dGDP, extent = None)
gr.get_poly_env()

mask = GDP_mask(dGDP)
gr.load_mask(mask)
t0 = time.time()
gr.load_data("GDP_PPP", first_dim = -1)
t1 = time.time()
print(f"{t1-t0:.2f}", "s.")

####################################################

def get_summed_value(index):
    poly = sh.gdf.iloc[index]["geometry"]
    output = gr.get_areafrac(poly,-1)
    subprocess.check_call(f"touch {d_log}{index}.log", shell = True)
    return output

####################################################
print("Start")

nbcore = multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=nbcore-1)
#

print("Number of polygons", sh.n_elements)

t0 = time.time()
P_list = pool.map(get_summed_value, np.arange(sh.n_elements))
t1 = time.time()
#
print(f"{t1-t0:.2f}", "s.")
#
out_gdf = sh.gdf.copy()
out_gdf["GDP"] = P_list
out_gdf.to_file(d_out_file, driver="GPKG")
#
t2 = time.time()
print(f"{t2-t1:.2f}", "s.")
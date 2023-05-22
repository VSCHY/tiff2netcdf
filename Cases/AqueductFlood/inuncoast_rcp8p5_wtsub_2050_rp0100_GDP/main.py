import sys
sys.path.append("../../../")
from src import Grid, GDP_mask
from src import Raster
import multiprocessing
import numpy as np
import xarray as xr
import tqdm

######################

daqflood = "/home/aschrapffer/Documents/tiff2netcdf/Cases/AqueductFlood/"
dGDP = "/home/aschrapffer/DATA/ECONOMIC/GDP_HDI_POP/GDP_PPP_1990_2015_5arcmin_v2.nc"
TEMP = "/home/aschrapffer/Documents/tiff2netcdf/TEMP/"

######################

name = "inuncoast_rcp8p5_wtsub_2050_rp0100_GDP"
dflood = "/home/aschrapffer/DATA/FLOOD/AqueductFloods/Coastal/inuncoast_rcp8p5_wtsub_2050_rp0100_0.tif"

######################

dire = daqflood + name + "/"
d_out_file = dire + "Output/" + name + ".nc"

######################

gr = Grid(dGDP, extent = None)

mask = GDP_mask(dGDP)
gr.load_mask(mask)
gr.select()

######################

rast = Raster(dflood, TEMP)

####################################################

def get_fraction_gt0(index):
    poly = gr.polylist[index]
    name = f"Riverine_100y_present_{index}"
    out = rast.get_fraction_threshold(poly, name, 0)
    return out 

def get_fraction_gt1(index):
    poly = gr.polylist[index]
    name = f"Riverine_100y_present_{index}"
    out = rast.get_fraction_threshold(poly, name, 1)
    return out 

####################################################

# RM LES LOG + recréer le dossier
block_size = 150000

n_part = gr.n_elements//block_size
nbcore = multiprocessing.cpu_count()

ds_dict = {}

####################################################
####################################################
# Greater than 0

print("Variable > 0m")
data_out_gt0 = gr.init_data_output()

for i in tqdm.tqdm(range(n_part+1)):
    pool = multiprocessing.Pool(processes=nbcore-1)
    if i<n_part:
        ind_land_indexes = np.arange(i*block_size, (i+1)*block_size)
    else:
        ind_land_indexes = np.arange(n_part*block_size, gr.n_elements)
    #
    P_list = pool.map(get_fraction_gt0, ind_land_indexes)
    #
    data_out_gt0 = gr.update_output(data_out_gt0, P_list, ind_land_indexes)
    pool.close()

data_out_gt0 = ma.masked_where(data_out_gt0 == -999., data_out_gt0)

da = gr.create_dataarray(data_out_gt0,  long_name = "Fraction with flood depth>0m", units = "-")

ds_dict["flooddepth_gt0"] = da

####################################################

print("Variable > 1m")
data_out_gt1 = gr.init_data_output()

for i in tqdm.tqdm(range(n_part+1)):
    pool = multiprocessing.Pool(processes=nbcore-1)
    if i<n_part:
        ind_land_indexes = np.arange(i*block_size, (i+1)*block_size)
    else:
        ind_land_indexes = np.arange(n_part*block_size, gr.n_elements)
    #
    P_list = pool.map(get_fraction_gt1, ind_land_indexes)
    #
    data_out_gt1 = gr.update_output(data_out_gt1, P_list, ind_land_indexes)
    pool.close()

data_out_gt1 = ma.masked_where(data_out_gt1 == -999., data_out_gt1)

da = gr.create_dataarray(data_out_gt1,  long_name = "Fraction with flood depth>1m", units = "-")

ds_dict["flooddepth_gt1"] = da

####################################################

ds = xr.Dataset(ds_dict,
                attrs = {})

ds.to_netcdf(d_out_file)

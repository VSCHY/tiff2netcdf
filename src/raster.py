import rasterio
from rasterio.mask import mask
import rioxarray
import numpy as np
import subprocess


class Raster:
    def __init__(self,d_raster, TEMP):
        self.d_raster = d_raster
        self.d_temp = TEMP


    def recut_poly(self, poly, name):
        with rasterio.open(self.d_raster) as src0:
            try:
                out_image, out_transform = mask(dataset=src0, shapes=[poly], crop=True)
                out_meta = src0.meta.copy()
            except:
                return 0
            

            out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})
            dname = self.d_temp + f"{name}.tif"

            with rasterio.open(dname, "w", **out_meta) as dest:
                dest.write(out_image)
            return 1

    def read_recut(self, name):
        dname = self.d_temp + f"{name}.tif"
        with rioxarray.open_rasterio(dname) as rds:
            lon = rds.x
            lat = rds.y
            data = rds.data[0,:,:]
            return lon, lat, data
        
    ##############################"
    
    def get_fraction(self, poly, name):
        res = self.recut_poly(poly, name)
        if res:
            _, _, data = self.read_recut(name)

            dfile = self.d_temp + f"{name}.tif"
            subprocess.check_call(f"rm {dfile}", shell = True)
            
            return np.sum(data)/data.size 
        else:
            return 0
    

            


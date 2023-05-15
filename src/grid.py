"""
Class to load GRID.
"""

from netCDF4 import Dataset
import numpy as np
from src import boxit
from shapely.geometry import Polygon
import xarray as xr
import numpy.ma as ma




class Grid:
   def __init__(self, d_grid, extent = None):
   
      nc = Dataset(d_grid, "r")
      
      var_list = list(nc.variables.keys())
      if "lon" in var_list:
          lon_name = "lon"
          lat_name = "lat"
      elif "longitude" in var_list:
          lon_name = "longitude"
          lat_name = "latitude"
          
      self.lon = nc.variables[lon_name][:]
      self.lat = nc.variables[lat_name][:]
      
      # Regular grid
      self.hdlon = np.mean(np.abs(np.diff(self.lon)))
      self.hdlat = np.mean(np.abs(np.diff(self.lat)))
      self.hd = 1* max(abs(self.hdlon), abs(self.hdlat))
      #
      self.lon = np.array([ll-360 if ll>180 else ll for ll in self.lon])

      self.extent = extent
      nc.close()      
   
   def reshape(self):
      if self.lon[0] >=0:
         i = np.where(self.lon<0)[0][0]
         self.lon = np.roll(self.lon, -i)
         ii = -i
      else:
         ii = 0  
      self.mask= np.roll(self.mask, ii, axis = 1)

      # Get the limits of the DOMAIN
      if self.extent is None:
          i0 = 0; j0 = 0
          j1 = self.lat.shape[0]
          i1 = self.lon.shape[0]
      else:
          i0 = np.argmin(np.abs(self.lon-(self.extent[0]-self.hd)))
          i1 = np.argmin(np.abs(self.lon-(self.extent[1]+self.hd)))
          i0, i1 = np.sort([i0,i1])
          j0 = np.argmin(np.abs(self.lat-(self.extent[2]-self.hd)))
          j1 = np.argmin(np.abs(self.lat-(self.extent[3]+self.hd)))
          j0, j1 = np.sort([j0,j1])
          print(i0,i1,j0,j1)
      self.limits = [j0, j1+1, i0, i1+1]

      self.lon = self.lon[i0:i1+1]
      self.lat = self.lat[j0:j1+1]
      self.nj = self.lat.shape[0]
      self.ni = self.lon.shape[0]
      #self.ind_land = [[j-j0,i-i0] for j in range(j1-j0) for i in range(i1-i0)]


   def load_mask(self, landseamask):
      self.mask = landseamask
      self.reshape()
      j0,j1,i0,i1 = self.limits
      self.mask = self.mask[j0:j1, i0:i1]
      
      A  = np.where(self.mask == 1)
      
      if len(A) == 1 :
         self.ind_land = [[j,i] for j in range(self.nj) for i in range(self.ni)]
      else:
         self.ind_land = [[A[0][i],A[1][i]] for i in range(len(A[0]))]
      self.get_corners()
      
   #
   def get_corners(self):
      self.centers = np.array([[self.lon[i], self.lat[j]] for j,i in self.ind_land])
   # 
   def get_poly(self, i):
      boxll = boxit(self.centers[i], self.hdlon, self.hdlat, 2)
      poly = Polygon(boxll)
      if not poly.is_valid:
         poly = poly.buffer(0)      
      return poly
   #
   def select(self):
      self.indices = np.arange(self.centers.shape[0])
      self.polylist = [self.get_poly(i) for i in self.indices]

      self.n_elements = self.indices.shape[0]

      #return indices, polylist, self.centers


   def create_output(self, P_list, varn, d_output):
      data_out = np.full((self.nj,self.ni), -999.)
      for indji, value in zip(self.ind_land,P_list):
         j,i = indji
         data_out[j,i] = value

      data_out = ma.masked_where(data_out == -999., data_out)
      print(data_out.dtype)

      lons = self.lon
      lats = self.lat


      xarlon = xr.IndexVariable("longitude", lons, attrs={"unit":"degrees_east","standard_name": "longitude", "long_name" : "longitude", "axis":"X"})
      xarlat = xr.IndexVariable("latitude", lats, attrs={"unit":"degrees_north","standard_name": "latitude", "long_name" : "latitude", "axis":"Y"})

      data_out = data_out.astype(np.float32)

      ds = xr.Dataset({
            varn: xr.DataArray(
                        data   = data_out,
                        dims   = ['latitude',"longitude"],
                        coords = {'latitude':xarlat,"longitude":xarlon},
                        attrs  = {"long_name":f"Fraction of cropland", 'units': '-'},
                        )},
                attrs = {})
      ds.to_netcdf(d_output)

"""
Contains a class to load GRID file (regular lonlat netCDF).
"""
from netCDF4 import Dataset
import numpy as np
from src import boxit
from shapely.geometry import Polygon
import xarray as xr
import numpy.ma as ma
from src import poly_env

##########################################################

class Grid:
   """
   Class to load netcdf file (regular grid).
   """
   def __init__(self, d_grid:str, extent = None):
      """
      Initialization of the class

      Inputs:
       - d_grid (str)  : location of the file
       - extent (list) : extent on which the file should be recut 
            (None by default i.e. full grid is loaded)
      """

      self.d_grid = d_grid
      nc = Dataset(d_grid, "r")
      self.ds = xr.open_dataset(self.d_grid)
      
      var_list = list(nc.variables.keys())
      if "lon" in var_list:
         self.lon_name = "lon"
         self.lat_name = "lat"
      elif "longitude" in var_list:
         self.lon_name = "longitude"
         self.lat_name = "latitude"
          
      self.lon = nc.variables[self.lon_name][:]
      self.lat = nc.variables[self.lat_name][:]

      if self.lat[1]<self.lat[0]:
         self.rev_lat = True
      else:
         self.rev_lat = False
      
      # Regular grid: regular lon/lat space
      self.hdlon = np.mean(np.abs(np.diff(self.lon)))
      self.hdlat = np.mean(np.abs(np.diff(self.lat)))
      self.hd = 1* max(abs(self.hdlon), abs(self.hdlat))

      # for issue with lon between 0 and 360
      self.lon = np.array([ll-360 if ll>180 else ll for ll in self.lon])
      self.extent = extent
      nc.close()    
      
      # 
      self.init_output()

   #############################  
   
   def reshape(self):
      """
      Reshape the data, to go from -180 to 180.
      """
      if self.lon[0] >=0:
         i = np.where(self.lon<0)[0][0]
         self.lon = np.roll(self.lon, -i)
         self.ii = -i
      else:
         self.ii = 0  
      self.mask= np.roll(self.mask, self.ii, axis = 1)

      # Get the limits of the DOMAIN (from extent)
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

      # Recut
      self.lon = self.lon[i0:i1+1]
      self.lat = self.lat[j0:j1+1]
      self.nj = self.lat.shape[0]
      self.ni = self.lon.shape[0]

   #############################
   
   def get_data(self, varn):
      data = self.ds[varn].values
      return data
   
   #############################

   def load_mask(self, landseamask):
      """
      Load the mask, numpy array filled with 0 and 1.
      The 1 are the area considered
      """
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
   
   #############################

   def load_data(self, varn:str, first_dim = None):
      """
      Load data, when the data from the netCDF is interpolated to shapefile.

      Inputs:
       - varn (str): name of the variable
       - first_dim (None or int): index of the data on first dimension when needed
         (default value is None if data is 2D).
      """
      self.varn = varn
      with Dataset(self.d_grid, "r") as nc:
         data = nc.variables[varn][:]
         data= np.roll(data, self.ii, axis = 1)

         j0,j1,i0,i1 = self.limits
         if first_dim is not None:
            self.data = data[first_dim, j0:j1, i0:i1]
         else:
            self.data = data[j0:j1, i0:i1]

   #############################

   def get_corners(self):
      """
      Get centers of the grid points.
      """
      self.centers = np.array([[self.lon[i], self.lat[j]] for j,i in self.ind_land])
   
   #############################

   def get_poly_index(self, i):
      """
      Get a polygon from the index of a grid point in the center list.
      """    
      return self.get_poly(self.centers[i])

   def get_poly(self, center):
      """
      Get a polygon from the index of a grid point in the center list.
      """
      boxll = boxit(center, self.hdlon, self.hdlat, 2)
      poly = Polygon(boxll)
      if not poly.is_valid:
         poly = poly.buffer(0)      
      return poly
   
   #############################
   
   def select(self):
      """
      Create list of corresponding polygons.
      """
      self.indices = np.arange(self.centers.shape[0])
      self.polylist = [self.get_poly_index(i) for i in self.indices]

      self.n_elements = self.indices.shape[0]
      #return indices, polylist, self.centers

   #############################

   def get_poly_env(self):
      """
      Load the poly_env module.
      """
      self.poly_env = poly_env(self.lon, self.lat)

   #############################

   def get_areafrac(self, shape, first_dim = None):
      """
      Get indexes of points and corresponding area fraction for a certain polygon.

      Input:
       - shape: polygon considered
      """
      
      minlon, minlat, maxlon, maxlat = shape.bounds
      lon0,lon1 = [float(minlon)-self.hd*0.6, float(maxlon)+self.hd*0.6]
      lat0,lat1 = [float(minlat)-self.hd*0.6, float(maxlat)+self.hd*0.6]

      if self.rev_lat:
         lat0,lat1 = (lat1,lat0)
      
      ds = self.ds.sel(longitude=slice(lon0,lon1), latitude=slice(lat0,lat1))
      lons = ds.coords[self.lon_name].values
      lats = ds.coords[self.lat_name].values
      if first_dim is not None:
         data = ds[self.varn].values[first_dim,:,:]
      else:
         data = ds[self.varn].values[:,:]
      NAN = np.isnan(data)      
      data = np.nan_to_num(data, nan = 0.0)
      nj,ni = data.shape

      indexes = [[j,i] for i in range(ni) for j in range(nj) if not NAN[j,i]]
      if len(indexes) == 0:
         del NAN
         return 0
      else: 
         centers = [np.array([lons[i], lats[j]]) for j,i in indexes]
         del lons; del lats
         polylist = [self.get_poly(cent) for cent in centers]   
         del centers   

         inter = np.array([self.poly_env.get_intersect(shape, point) for point in polylist])
         del polylist
         inter = inter[:,1]
         data = [data[j,i] for j,i in indexes]
         del indexes
         output = [frac*d for frac, d in zip(inter,data)]
         del data
         ds.close()
         return np.sum(output)
   
   #############################

   def init_output(self):

      self.xarlon = xr.IndexVariable("longitude", self.lon, attrs={"unit":"degrees_east","standard_name": "longitude", "long_name" : "longitude", "axis":"X"})
      self.xarlat = xr.IndexVariable("latitude", self.lat, attrs={"unit":"degrees_north","standard_name": "latitude", "long_name" : "latitude", "axis":"Y"})   

   def create_output(self, P_list, varn:str, d_output:str, ind_land = None):
      """
      Create netCDF output.

      Input:
       - P_list: list of values
       - varn: desired variable name
       - d_output: directory and name of the output file
      """
      data_out = np.full((self.nj,self.ni), -999.)
      if ind_land is None:
          ind_land = self.ind_land
          
      for indji, value in zip(ind_land, P_list):
         j,i = indji
         data_out[j,i] = value
      data_out = ma.masked_where(data_out == -999., data_out)
      data_out = data_out.astype(np.float32)
      
      ds = xr.Dataset({
            varn: xr.DataArray(
                        data   = data_out,
                        dims   = ['latitude',"longitude"],
                        coords = {'latitude':self.xarlat,"longitude":self.xarlon},
                        attrs  = {"long_name":f"Fraction of cropland", 'units': '-'},
                        )},
                attrs = {})
      ds.to_netcdf(d_output)
      ds.close()
      del ds; del data_out

   #############################

   def init_data_output(self, fill_number = -999.):
      data_out = np.full((self.nj,self.ni), fill_number)
      return data_out

   def update_output(self, data_out, P_list, ind_land_indexes = None):
      if ind_land_indexes is None:
         ind_land_indexes = np.arange(len(self.ind_land))
          
      for indji, value in zip(ind_land_indexes, P_list):
         j,i = self.ind_land[indji]
         data_out[j,i] = value
      return data_out

   def create_dataarray(self, data_out, long_name = "", units = "-"):
      da = xr.DataArray(
            data   = data_out,
            dims   = ['latitude',"longitude"],
            coords = {'latitude':self.xarlat,"longitude":self.xarlon},
            attrs  = {"long_name":long_name, 'units': units},
            )
      return da
   
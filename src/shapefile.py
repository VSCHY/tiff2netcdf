import geopandas as gpd
from shapely.geometry import Point

class Shapefile:
    def __init__(self, d_file):
        self.gdf = gpd.read_file(d_file)
        self.n_elements = self.gdf.shape[0]
        self.bounds = self.gdf.bounds
        
    def get_poly_from_lonlat(self, lon0, lat0):
        bb = self.bounds.copy()
        bb = bb[lon0 < bb["maxx"]]
        bb = bb[lon0 > bb["minx"]]
        bb = bb[lat0 < bb["maxy"]]
        bb = bb[lat0 > bb["miny"]]
        
        gdf_out = self.gdf.loc[bb.index]
        if gdf_out.shape[0] == 1:
            return gdf_out
        else:
            if gdf_out.shape[0] > 1:
                gdf_out = gdf_out[gdf_out.contains(Point(lon0,lat0))] 
            if gdf_out.shape[0] == 0:
                return None
            else:
                return gdf_out
            
    
  

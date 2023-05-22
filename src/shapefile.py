import geopandas as gpd


class Shapefile:
    def __init__(self, d_file):
        self.gdf = gpd.read_file(d_file)
        self.n_elements = self.gdf.shape[0]
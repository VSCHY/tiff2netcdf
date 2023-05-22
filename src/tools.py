
from netCDF4 import Dataset
from numba import jit
import numpy as np
from shapely.geometry import Polygon

@jit(nopython = True)
def boxit(cent, dlon, dlat, polyres) :
    boxll=[]
    loninc = dlon/polyres
    latinc = dlat/polyres
    # Side 1
    for pts in range(polyres) :
        boxll.append([cent[0]-dlon/2.0+loninc*pts, cent[1]+dlat/2.0])
    # Side 2
    for pts in range(polyres) :
        boxll.append([cent[0]+dlon/2.0, cent[1]+dlat/2.0-latinc*pts])
    # Side 3
    for pts in range(polyres) :
        boxll.append([cent[0]+dlon/2.0-loninc*pts, cent[1]-dlat/2.0])
    # Side 4
    for pts in range(polyres) :
        boxll.append([cent[0]-dlon/2.0, cent[1]-dlat/2.0+latinc*pts])
    # Close
    boxll.append(boxll[0])
    return boxll


def boxit_lonslats(cent, dlon, dlat, polyres):
    lons = np.zeros(polyres*4, dtype = np.float32)
    lats = np.zeros(polyres*4, dtype = np.float32)
    
    loninc = dlon/polyres
    latinc = dlat/polyres
    
    i = 0
    # Side 1
    for pts in range(polyres) :
        lats[i] = cent[1]+dlat/2.0
        lons[i] = cent[0]-dlon/2.0+loninc*pts
        i+=1
    # Side 2
    for pts in range(polyres) :
        lats[i] = cent[1]+dlat/2.0-latinc*pts
        lons[i] = cent[0]+dlon/2.0
        i+=1
    # Side 3
    for pts in range(polyres) :
        lats[i] = cent[1]-dlat/2.0
        lons[i] = cent[0]+dlon/2.0-loninc*pts
        i+=1
    # Side 4
    for pts in range(polyres) :
        lats[i] = cent[1]-dlat/2.0+latinc*pts
        lons[i] = cent[0]-dlon/2.0
        i+=1

    return lons, lats

def ERA5_mask(dlandsea):
    with Dataset(dlandsea, "r") as nc:
        mask  = nc["lsm"][0,:,:]
        mask[mask>0] = 1
    return mask

def GDP_mask(dGDP):
    with Dataset(dGDP, "r") as nc:
        mask  = nc["GDP_PPP"][0,:,:]
        mask[mask!=-9] = 1
        mask[mask==-9] = 0
    return mask

def POP_mask(dPOP):
    with Dataset(dPOP, "r") as nc:
        mask  = nc["UN WPP-Adjusted Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes"][0,:,:]
        mask = mask.mask.astype(np.int32)
        mask = 1-mask
    return mask

class poly_env:
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat
        self.dlon = np.unique(np.diff(lon))
        self.dlon = self.dlon[self.dlon!=0]
        self.dlon = np.abs(self.dlon[0])
        self.dlat = np.unique(np.diff(lat))
        self.dlat = self.dlat[self.dlat!=0]
        self.dlat = np.abs(self.dlat[0])

    def get_intersect(self, shape, p):
        if shape.intersects(p):
            af = shape.intersection(p).area / p.area
            return [True, af]
        else:
            return [False, 0]
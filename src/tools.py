
from netCDF4 import Dataset
from numba import jit

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


def ERA5_mask(dlandsea):
    with Dataset(dlandsea, "r") as nc:
        mask  = nc["lsm"][0,:,:]
        mask[mask>0] = 1
    return mask

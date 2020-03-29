#################################################################################################
# Sample script to open the high resolution atmospheric moisture tracking data based on UTrack  #
#                                                                                               #
# The dataset (Tuinenburg et al., 2020) is available from the PANGAEA archive in both the full  #
# resolution version at 0.5° resolution and in a lower resolution version at 1.0° resolution    #
# (https://doi.pangaea.de/10.1594/PANGAEA.912710).                                              #
# Contact: Obbe Tuinenburg, O.A.Tuinenburg <at> uu.nl                                                #
#################################################################################################

# Note that the following script reads a set of footprints for one location.
# If you need more than one location, it is much more efficient to read all the data from the netcdf file
# at once, instead of one by one.

from scipy.io.netcdf import *
from netCDF4 import Dataset
from numpy import *
import matplotlib.pyplot as plt
import sys

# Please provide the month number on the command line, otherwise take July
try:
    month=int(sys.argv[1])
except:
    month=7
f=Dataset('utrack_climatology_0.5_'+str(month).zfill(2)+'.nc')
etfile=Dataset('ERA5_ET_0.5.nc') # Netcdf file with monthly ERA5 ET
ET=-1000.0*mean(etfile.variables['ET'][arange(month-1,etfile.variables['ET'][:].shape[0],12)],0)

lats=arange(90,-90,-0.5)
lons=arange(0,360,0.5)

def get_closest_index(lats,lat):
        import operator
        lat_index, min_value = min(enumerate(abs(lats-lat)), key=operator.itemgetter(1))
        return lat_index

def get_footprints(latitude,longitude):
    # Determine the closest indices for the location
    latidx=get_closest_index(lats,latitude)
    lonidx=get_closest_index(lons,longitude)

    # Determine the forward footprint, where the ET from the input location will subsequently rain out.
    fp=f.variables['moisture_flow'][latidx,lonidx]
    fp=fp*-0.1
    fp=e**fp
    fp[fp==1]=0
    forward_fp=fp/sum(fp)

    # Determine the backward footprint, where the precipitation falling in the input location has previously evaporated
    fp=f.variables['moisture_flow'][:,:,latidx,lonidx]
    fp=fp*-0.1
    fp=e**fp
    fp[fp==1]=0
    fp=ET*fp
    backward_fp=fp/sum(fp)
    return forward_fp,backward_fp

# Footprints in Utrecht
forward_footprint,backward_footprint=get_footprints(52,5.1)

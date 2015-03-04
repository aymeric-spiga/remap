#! /usr/bin/env python
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import ppplot

### READ NETCDF
f = NetCDFFile('yorgl.nc')
lons = f.variables['lon'][:]
lats = f.variables['lat'][:]
z = f.variables['ps'][:]

#### FIND UNIQUE LAT/LON
#zelat,colat=np.unique(np.round(lats,decimals=8),return_inverse=True)
#zelon,colon=np.unique(np.round(lons,decimals=8),return_inverse=True)
#### RESHAPE TO DIMENSIONS SUGGESTED BY UNIQUE LAT/LON
#shp = (zelat.shape[0],zelon.shape[0])
##print lats.shape,shp[0]*shp[1]

### GUESS SIZE
N = np.int(np.sqrt(lons.shape[0]/2))
shp = (N,2*N)

### RESHAPE
lat = np.reshape(lats,shp)
lon = np.reshape(lons,shp)
ff = np.reshape(z,shp)

print "plot"

### PLOT
import ppplot
pl = ppplot.plot2d()
pl.f = ff
pl.x = lon
pl.y = lat
pl.proj = "ortho"
pl.blat = 45.
pl.makeshow()
exit()

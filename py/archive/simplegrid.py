#! /usr/bin/env python

import netCDF4 as nc
import sys
import math
import numpy as np


def from_reduced(N,M):
  #"N elements from south to north and N elements around equator "
  cells_lon = []
  cells_lat = []
  for i in range(-M/2,M/2):
    lat1 = 180.0/M*i
    lat2 = 180.0/M*(i+1)
    #print "yorgl",i,lat1,lat2
    for j in range(-N/2,N/2):
      lon1 = 360.0*j/N
      lon2 = 360.0*(j+1)/N
      if (i == -M/2): 
        print "yorgl",j,lon1,lon2
      bounds_lon = [lon1, lon1, lon2, lon2]
      bounds_lat = [lat1, lat2, lat2, lat1]
      bounds_lon.append(bounds_lon[0])
      bounds_lat.append(bounds_lat[0])
      cells_lon.append(bounds_lon) 
      cells_lat.append(bounds_lat)
  return np.array(cells_lon), np.array(cells_lat)




#for N in [64, 128, 256, 512]:
for N in [32,256]:
	filename = "simple" + str(N) + ".nc"

	print "Generating: N =", N 
	#lon, lat = from_reduced(N*2,N)
        lon, lat = from_reduced(N,N)

	print lon.shape[0], "cells -> writing as ", filename

	f = nc.Dataset(filename,'w')
	f.createDimension('n_vert', 5)
	f.createDimension('n_cell', lon.shape[0])

	var = f.createVariable('lat', 'd', ('n_cell'))
	var.setncattr("long_name", "latitude")
	var.setncattr("units", "degrees_north")
	var.setncattr("bounds", "bounds_lat")
	var[:] = np.zeros(lon.shape[0])
	var = f.createVariable('lon', 'd', ('n_cell'))
	var.setncattr("long_name", "longitude")
	var.setncattr("units", "degrees_east")
	var.setncattr("bounds", "bounds_lon")
	var[:] = np.zeros(lon.shape[0])

	var = f.createVariable('bounds_lon', 'd', ('n_cell','n_vert'))
	var[:] = lon
	var = f.createVariable('bounds_lat', 'd', ('n_cell','n_vert'))
	var[:] = lat
	var = f.createVariable('val', 'd', ('n_cell'))
	var.setncattr("coordinates", "lon lat")
	var[:] = np.arange(lon.shape[0])
	f.close()


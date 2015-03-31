#! /usr/bin/env python

import time
import numpy as np

def ll(N,M=None):
  # MAKE A SIMPLE 2NxN LON/LAT GRID 
  # use array operations for speed-up
  # ---- N MUST BE AN EVEN NUMBER
  # Author: A. Spiga
  stime = time.time()
  N = int(N) 
  if M is None: M = 2*N
  deg = 180.0 #whole sphere
  #deg = 90.0  #works
  #deg = 60.0  #works
  #deg = 45.0  #does not work
  #deg = 30.0  #does not work
  lat1 = (deg/N)*np.arange(-N/2,N/2)
  lat2 = (deg/N)*np.arange(1-N/2,1+N/2)
  lon1 = (2.*deg/M)*np.arange(-M/2,M/2)
  lon2 = (2.*deg/M)*np.arange(1-M/2,1+M/2)
  ## make 2D versions of coordinates
  tablon1,tablat1 = np.meshgrid(lon1,lat1)
  tablon2,tablat2 = np.meshgrid(lon2,lat2)
  ## make tabs 1D only -- to prepare cells
  tablon1 = np.ravel(tablon1) ; tablon2 = np.ravel(tablon2)
  tablat1 = np.ravel(tablat1) ; tablat2 = np.ravel(tablat2)
  ## stack tabs together to indicate bounds and center for each point
  tabbounds_lon = np.dstack((tablon1,tablon1,tablon2,tablon2,tablon1))
  tabbounds_lat = np.dstack((tablat1,tablat2,tablat2,tablat1,tablat1))
  ## delete dimension of rank 1 in previous operation
  cells_lon = np.squeeze(tabbounds_lon)
  cells_lat = np.squeeze(tabbounds_lat)
  ## output bounds lon, bounds lat, center lon, center lat
  return cells_lon, cells_lat, tablon1, tablat1

def apply_weights(src_val_glo,A):
  src_val = np.hstack(src_val_glo)
  dst_val = A*src_val
  ## trick to remove multiple cells
  ## ponderate by number of multiple cells
  ## this number is determined 
  ## by applying A to a matrix full of 1s
  src_val = src_val*0. + 1.
  dst_val = dst_val / (A*src_val)
  return dst_val

def compute_distribution(ncell,rank,size):
  "Returns the local number and starting position in global array."
  if rank < ncell % size:
    return ncell//size + 1, \
           (ncell//size + 1)*rank
  else:
    return ncell//size, \
           (ncell//size + 1)*(ncell%size) + (ncell//size)*(rank - ncell%size)


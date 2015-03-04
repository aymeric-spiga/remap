#! /usr/bin/env python

##
## REMAP
##
## rewrite + generalization + options + optimization: A. Spiga Feb 2015
##

import netCDF4 as nc
import ctypes as ct
import numpy as np
import os
import time
from optparse import OptionParser ### TBR by argparse
from remap_func \
  import ll,apply_weights,compute_distribution


## parallel or not
parallel=False
parallel=True
if parallel: from mpi4py import MPI

##########
## DICT ##
##########
grid_types = {
"dynamico:mesh": {"lon_name": "bounds_lon_i","lat_name": "bounds_lat_i","pole": [0,0,0]},
"dynamico:vort": {"lon_name": "bounds_lon_v","lat_name": "bounds_lat_v","pole": [0,0,0]},
"dynamico:restart": {"lon_name": "lon_i_vertices","lat_name": "lat_i_vertices","pole": [0,0,0]},
"test:polygon": {"lon_name": "bounds_lon","lat_name": "bounds_lat","pole": [0,0,0]},
"test:latlon": {"lon_name": "bounds_lon","lat_name": "bounds_lat","pole": [0,0,1]},
"ll": {"func": ll,"pole": [0,0,1]},
}
interp_types = {"FV1": 1,"FV2": 2}

######################################
# define parser with version and usage 
######################################
parser = OptionParser()
parser.version = \
'''
REMAPPER 
`remap.py -h` for usage and options
'''
##
parser.usage = \
'''
remap.py [options] srcfile dstfile
--> srcfile is an input netCDF file
--> dstfile is either a number for simple latlon grid
                   or a netCDF file indicating destination grid
'''
##
parser.add_option('-W','--onlyweights',action='store_true',dest='onlyweights',default=False,\
  help="only compute weights")
parser.add_option('-F','--forceweights',action='store_true',dest='forceweights',default=False,\
  help="force computing of weights")
parser.add_option('-S','--srctype',action='store',dest='srctype',default="test:polygon",\
  help="grid type of source")
parser.add_option('-D','--dsttype',action='store',dest='dsttype',default="ll",\
  help="grid type of destination")
parser.add_option('-o','--outfile',action='store',dest='outfile',type="string",default="outremap.nc",\
  help="output file")
parser.add_option('-R','--reshaped',action='store_true',dest='reshaped',default=False,\
  help="output reshaped fields on a 2D grid")
parser.add_option('-i','--interp',action='store',dest='interp',type="string",default="FV1",\
  help="interpolation method (FV1 FV2) first order conservative Finite Volume")
parser.add_option('-v','--var2d',action='append',dest='var2d',type="string",default=None,\
  help="2D field [append is possible]")
parser.add_option('-V','--var3d',action='append',dest='fieldchar',type="string",default=None,\
  help="3D field [append is possible]")
parser.add_option('-Z','--vert',action='store',dest='vertchar',type="string",default="presnivs",\
  help="vertical coordinate")
parser.add_option('-z','--level',action='append',dest='z',type="int",default=None,\
  help="choose vertical indexes to be interpolated [append is possible]")
parser.add_option('-P','--plot',action='store_true',dest='plot',default=False,\
  help="plot fields")
##
(opt,args) = parser.parse_args()
if (len(args) == 0): parser.print_version() ; exit()
## SRCFILE
if len(args) != 2: parser.print_usage() ; exit(2)
else: srcfile = args[0] ; dstfile = args[1]
## SRCTYPE
try:
    srctype = grid_types[opt.srctype]
except KeyError:
    print "Error: srctype needs to be one of the following: " + " ".join(grid_types.keys()) + "."
    exit(2)
## DSTTYPE
try:
    dsttype = grid_types[opt.dsttype]
except KeyError:
    print "Error: dsttype needs to be one of the following: " + " ".join(grid_types.keys()) + "."
    exit(2)
## FIELDCHAR
if opt.fieldchar is None: opt.fieldchar = ["temp"]
if opt.var2d is None: opt.var2d = ["ps"]
## NO SPECIFIC OPERATION NEEDED
vertchar = opt.vertchar
interp = opt.interp

##
## test if we have to compute weights
##
wfile = srcfile+"_weights"
if "func" in dsttype: wfile = wfile + "_" + dstfile
wfile = wfile+'.nc'
if opt.forceweights: computeweights = False
else: computeweights = not(os.path.isfile(wfile))

##
## test if we have to compute barycentres
##
if "func" in dsttype: computebary = False
else: computebary = True

####
#### LOAD remap LIBRARY
####
remap = ct.cdll.LoadLibrary(os.path.realpath('libmapper.so'))

###
### MAIN PROGRAM
###

remap.mpi_init()
rank = remap.mpi_rank()
size = remap.mpi_size()
#print rank, "/", size

############
### GRID ###
############
print "**** GRID ****"
print "Get grids either from files or from computations"
stime = time.time()

if "reader" in srctype:
	src_lon, src_lat = srctype["reader"](srcfile)
else:
	src = nc.Dataset(srcfile)
	# the following two lines do not perform the actual read
	# the file is read later when assigning to the ctypes array
	# -> no unnecessary array copying in memory
	src_lon = src.variables[srctype["lon_name"]]
	src_lat = src.variables[srctype["lat_name"]]

if "reader" in dsttype:
	dst_lon, dst_lat = dsttype["reader"](dstfile)
elif "func" in dsttype:
        dst_lon, dst_lat, dst_centre_lon, dst_centre_lat = dsttype["func"](dstfile)
else:
	dst = nc.Dataset(dstfile)
	dst_lon = dst.variables[dsttype["lon_name"]]
	dst_lat = dst.variables[dsttype["lat_name"]]

## prepare dimensions and arrays for later use in library and computations
dst_ncell, dst_nvert = dst_lon.shape
dst_ncell_loc, dst_loc_start = compute_distribution(dst_ncell,rank,size)
dstpole = (ct.c_double * (3))() ; dstpole[:] = dsttype["pole"]
c_dst_ncell = ct.c_int(dst_ncell_loc)
c_dst_nvert = ct.c_int(dst_nvert)
order = ct.c_int(interp_types[interp])
if computeweights or computebary:
  print "convert and reshape to C-type arrays for remap library"
  ## lon
  c_dst_lon = (ct.c_double * (dst_ncell_loc*dst_nvert))()
  zelen = len(c_dst_lon)
  c_dst_lon[:] = nc.numpy.reshape(dst_lon[dst_loc_start:dst_loc_start+dst_ncell_loc,:], (zelen,1))
  ## lat
  c_dst_lat = (ct.c_double * (dst_ncell_loc*dst_nvert))()
  c_dst_lat[:] = nc.numpy.reshape(dst_lat[dst_loc_start:dst_loc_start+dst_ncell_loc,:], (zelen,1))

gtime = time.time() - stime

###############
### WEIGHTS ###
###############
print "**** WEIGHTS ****"
stime = time.time()
### -- if weight file does not exist, calculate weights and create file
### -- if weight file does exist, read weights
if rank == 0:
 if computeweights:
    print "Calling remap library to compute weights."
    c_nweight = ct.c_int()
    ## convert to C-type arrays for remap library
    src_ncell, src_nvert = src_lon.shape
    src_ncell_loc, src_loc_start = compute_distribution(src_ncell,rank,size)
    c_src_lon = (ct.c_double * (src_ncell_loc*src_nvert))()
    zelen = len(c_src_lon)
    c_src_lon[:] = nc.numpy.reshape(src_lon[src_loc_start:src_loc_start+src_ncell_loc,:], (zelen,1))
    c_src_lat = (ct.c_double * (src_ncell_loc*src_nvert))()
    c_src_lat[:] = nc.numpy.reshape(src_lat[src_loc_start:src_loc_start+src_ncell_loc,:], (zelen,1))
    srcpole = (ct.c_double * (3))() ; srcpole[:] = srctype["pole"]
    c_src_ncell = ct.c_int(src_ncell_loc)
    c_src_nvert = ct.c_int(src_nvert)
    ##
    print "remap_get_num_weights"
    remap.remap_get_num_weights(c_src_lon, c_src_lat, c_src_nvert, c_src_ncell, srcpole,
               c_dst_lon, c_dst_lat, c_dst_nvert, c_dst_ncell, dstpole,
               order, ct.byref(c_nweight))
    ##
    nwgt = c_nweight.value
    c_weights = (ct.c_double * nwgt)()
    c_dst_idx = (ct.c_int * nwgt)()
    c_src_idx = (ct.c_int * nwgt)()
    ##
    print "remap_get_weights"
    remap.remap_get_weights(c_weights, c_src_idx, c_dst_idx)
    ##
    if parallel:
      wgt_glo     = MPI.COMM_WORLD.gather(c_weights[:])
      src_idx_glo = MPI.COMM_WORLD.gather(c_src_idx[:])
      dst_idx_glo = MPI.COMM_WORLD.gather(c_dst_idx[:])
    else:
      wgt_glo     = c_weights[:]
      src_idx_glo = c_src_idx[:]
      dst_idx_glo = c_dst_idx[:]
    ### change lists to numpy arrays to be saved and used in calculations
    wgt_glo     = np.hstack(wgt_glo)
    src_idx_glo = np.hstack(src_idx_glo)
    dst_idx_glo = np.hstack(dst_idx_glo)
    ### create netCDF file
    nwgt_glo = wgt_glo.size
    print "Writing", nwgt_glo, "weights to netCDF-file '" + wfile + "'."
    f = nc.Dataset(wfile,'w')
    f.createDimension('n_src',    src_ncell)
    f.createDimension('n_dst',    dst_ncell)
    f.createDimension('n_weight', nwgt_glo)
    var = f.createVariable('src_idx', 'i', ('n_weight')) ; var[:] = src_idx_glo
    var = f.createVariable('dst_idx', 'i', ('n_weight')) ; var[:] = dst_idx_glo
    var = f.createVariable('weight',  'd', ('n_weight')) ; var[:] = wgt_glo
    f.close()
 else:
    print "Reading weights from netCDF file "+wfile
    f = nc.Dataset(wfile)
    src_idx_glo = f.variables['src_idx'][:]
    dst_idx_glo = f.variables['dst_idx'][:]
    wgt_glo = f.variables['weight'][:]
    f.close()
wtime = time.time() - stime

#############
### REMAP ###
#############
stime = time.time()
if not opt.onlyweights:

    print "**** REMAP ****"

    ### Barycentres and areas if needed
    if computebary:
      print 'Get barycentres and areas'
      ##
      c_centre_lon = (ct.c_double * dst_ncell_loc)()
      c_centre_lat = (ct.c_double * dst_ncell_loc)()
      c_areas      = (ct.c_double * dst_ncell_loc)()
      remap.remap_get_barycentres_and_areas(c_dst_lon, c_dst_lat, c_dst_nvert, c_dst_ncell, dstpole,
		c_centre_lon, c_centre_lat, c_areas)
      ##
      if parallel:
        dst_areas_glo = MPI.COMM_WORLD.gather(np.array(c_areas[:]))
        dst_centre_lon_glo = MPI.COMM_WORLD.gather(np.array(c_centre_lon[:]))
        dst_centre_lat_glo = MPI.COMM_WORLD.gather(np.array(c_centre_lat[:]))
      else:
        dst_areas_glo = np.array(c_areas[:])
        dst_centre_lon_glo = np.array(c_centre_lon[:])
        dst_centre_lat_glo = np.array(c_centre_lat[:])
      ##
      if rank == 0:
          dst_areas = np.hstack(dst_areas_glo)
          dst_centre_lon = np.hstack(dst_centre_lon_glo)
          dst_centre_lat = np.hstack(dst_centre_lat_glo)

    ### determine vertical levels
    presnivs=src.variables[vertchar]
    nz = len(presnivs)
    if opt.z is None:
      vertrange = range(nz)
    else:
      vertrange = opt.z ; nz = len(vertrange)

    ### Prepare netCDF file for write
    f = nc.Dataset(opt.outfile,'w')
    ### first treat vertical coordinates
    f.createDimension(vertchar, nz)
    var = f.createVariable(vertchar, 'd', (vertchar))
    var.setncattr("long_name", "vertical coordinate")
    var.setncattr('axis', 'Z')
    if opt.z is None: var[:] = presnivs[:]
    else: var[:] = presnivs[opt.z]
    ### second, horizontal coordinates
    ### -- two modes: based on cells (default), or based on lat/lon (if reshaped=True)
    if opt.reshaped:
        ps = np.zeros( dst_ncell )
        temp = np.zeros( (nz,dst_ncell) )
        ##
        N = np.int(np.sqrt(dst_ncell/2))
        nt = 1
        shp = (nt,N,N*2)
        shp3 = (nt,nz,N,N*2)
        shphor = (N,N*2)
        ##
        f.createDimension('lon', N*2)
        f.createDimension('lat', N)
        f.createDimension('Times', nt)
        ##
        var = f.createVariable('lon', 'd', ('lon'))
        var.setncattr("long_name", "longitude")
        var.setncattr("units", "deg north")
        var[:] = np.unique(dst_centre_lon)[:]
        ##
        var = f.createVariable('lat', 'd', ('lat'))
        var.setncattr("long_name", "latitude")
        var.setncattr("units", "deg east")
        var[:] = np.unique(dst_centre_lat)[:]
        ##
        var = f.createVariable('Times', 'd', ('Times'))
        var[:] = [0.]
    else:
#	nq=src.dimensions['nq']
	f.createDimension('nvert', dst_nvert)
	f.createDimension('cell', dst_ncell)
#	f.createDimension('nq', len(nq))

	var = f.createVariable('lat', 'd', ('cell'))
	var.setncattr("long_name", "latitude")
	var.setncattr("units", "degrees_north")
	var.setncattr("bounds", "bounds_lat")
	var[:] = dst_centre_lat
	var = f.createVariable('lon', 'd', ('cell'))
	var.setncattr("long_name", "longitude")
	var.setncattr("units", "degrees_east")
	var.setncattr("bounds", "bounds_lon")
	var[:] = dst_centre_lon

	var = f.createVariable('bounds_lon', 'd', ('cell','nvert'))
	var[:] = dst_lon
	var = f.createVariable('bounds_lat', 'd', ('cell','nvert'))
	var[:] = dst_lat

    ########################################################################
    ### THE A MATRIX TO CHANGE COORDINATES ###
    if rank == 0:
        from scipy import sparse
        A = sparse.csr_matrix(sparse.coo_matrix((wgt_glo,(dst_idx_glo,src_idx_glo))))
    ########################################################################

    # 2D FIELD
    for var2d in opt.var2d:
      print "remapping...",var2d
      src_val_loc = src.variables[var2d]
      dim = len(src_val_loc.shape)
      if dim == 2: tab_loc = np.array(src_val_loc[0,:])
      elif dim == 1: tab_loc = np.array(src_val_loc[:])
      if parallel:
        src_val_glo = MPI.COMM_WORLD.gather(tab_loc)
      else:
        src_val_glo = tab_loc
      if rank == 0:
          dst_val = apply_weights(src_val_glo,A)
          if not opt.reshaped:
            ps = f.createVariable(var2d, 'd', ('cell'))
            ps.setncattr("coordinates", "lon lat")
          ps[:] = dst_val
      if opt.reshaped:
        print "reshaping and writing...",var2d
        var = f.createVariable(var2d, 'd', ('Times','lat','lon'))
        var[:,:,:] = np.reshape(ps,shp)

    # 3D FIELD
    for fieldchar in opt.fieldchar:
      print "remapping... %s with %i levels" % (fieldchar,nz)
      src_val_loc = src.variables[fieldchar]
      dim = len(src_val_loc.shape)
      if not opt.reshaped:
        temp = f.createVariable(fieldchar, 'd', (vertchar,'cell'))
        temp.setncattr("coordinates", "presnivs lon lat")
      countlev=0
      tmptime = time.time()
      for l in vertrange:
        if dim == 3: tab_loc = np.array(src_val_loc[0,l,:])
        elif dim == 2: tab_loc = np.array(src_val_loc[l,:])
        if parallel:
          src_val_glo = MPI.COMM_WORLD.gather(tab_loc)
        else:
          src_val_glo = tab_loc
        if rank == 0:
            dst_val = apply_weights(src_val_glo,A)
            temp[countlev,:] = dst_val
        test = time.time() - tmptime
        if test > 5.: 
          print "5s elapsed. done up to vertical level %i" % (l)
          tmptime = time.time()
        countlev = countlev+1
      if opt.reshaped:
        print "reshaping and writing...",fieldchar
        var = f.createVariable(fieldchar, 'd', ('Times',vertchar,'lat','lon'))
        var[:,:,:,:] = np.reshape(temp,shp3)
        print "...done"
 
rtime = time.time() - stime
f.close()

############
### PLOT ###
############
if not opt.onlyweights:
 if ("func" in dsttype) and (opt.plot):
  import ppplot
  print "**** PLOT ****"
  ### GUESS SIZE
  N = np.int(np.sqrt(dst_ncell/2))
  shp = (N,N*2)
  shp3 = (nz,N,N*2)
  ### PLOT
  pl = ppplot.plot2d()
  ### RESHAPE (should be at no cost) and ASSIGN to PLOT
  field = ps
  field = temp[0,:]
  pl.f = np.reshape(field,shp)
  pl.y = np.reshape(dst_centre_lat,shp)
  pl.x = np.reshape(dst_centre_lon,shp)
  ### PLOT SETTINGS and MAKE
  pl.proj = "cyl"#"ortho"
  pl.blat = 45.
  pl.makeshow()
  ### A ZONAL SECTION
  tab = np.reshape(temp,shp3)
  pl.f = np.mean(tab,axis=2)
  pl.x = None ; pl.y = None
  pl.makeshow()

###########
### END ###
###########
if not "reader" in srctype:
	src.close()
if not "reader" in dsttype:
  if not "func" in dsttype:
	dst.close()

print "**** TIMES ****"
print "GRID %.2f sec // WEIGHTS %.2f sec // REMAP %.2f sec // TOTAL %.2f sec" % (gtime,wtime,rtime,gtime+wtime+rtime)

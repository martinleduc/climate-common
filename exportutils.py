from IPython import embed
import numpy as np
# the Scientific Python netCDF 3 interface
# http://dirac.cnrs-orleans.fr/ScientificPython/
#from Scientific.IO.NetCDF import NetCDFFile as Dataset
# the 'classic' version of the netCDF4 python interface
# http://code.google.com/p/netcdf4-python/
from netCDF4 import Dataset
from numpy import arange, dtype # array module from http://numpy.scipy.org

def mat2esriascii(lon,lat,data,out):
    '''Export a 2D matrix to ESRI ASCII raster format for GIS.'''

    data=np.flipud(data)

    fi=open(out,'w')

    ncols=data.shape[1]
    nrows=data.shape[0]
    xllcenter=lon[0]
    yllcenter=lat[0]
    cellsize=lon[1]-lon[0]
    nodata_value= -9999999999999999


    fi.write('ncols '+str(ncols)+'\n')
    fi.write('nrows '+str(nrows)+'\n')
    fi.write('xllcenter '+str(xllcenter)+'\n')
    fi.write('yllcenter '+str(yllcenter)+'\n')
    fi.write('cellsize '+str(cellsize)+'\n')
    fi.write('nodata_value -9999999999'+'\n')

    for jj in range(data.shape[0]):
        fi.write(' '.join(str(data[jj,:]).strip('[]').split(','))+'\r\n')
    
    fi.close()
    print '.asc file written successfully.'


def mat2asciicolumns(lon,lat,data,field,out,sep='     ;     '):
    '''Export a 2D matrix to ESRI ASCII raster format for GIS.'''

    nlon=len(lon)
    nlat=len(lat)

    fi=open(out,'w')

    line="%14s %s %14s %s %14s \r\n" % ('longitude',sep,'latitude',sep,field)
    fi.write(line)
    fi.write('\r\n')
    
    for jj in range(nlat):
        for ii in range(nlon):
            line="%14f %s %14f %s %14f \r\n" % (lon[ii],sep,lat[jj],sep,data[jj,ii])
            fi.write(line)



            
    fi.close()
    print '.asc file written successfully.'

    

def mat2netCDF4(lon_in,lat_in,data_in,varname,ncpath,attrib=''):
    ncfile = Dataset(ncpath,'w')

    ny,nx=data_in.shape

    # create dimensions.
    ncfile.createDimension('x',nx)
    ncfile.createDimension('y',ny)
    #ncfile.createDimension('t',nt)

    # create the variables
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    # data = ncfile.createVariable('tas',dtype('float32').char,('x','y'))
    # lon = ncfile.createVariable('lon',dtype('double').char,('x'))
    # lat = ncfile.createVariable('lat',dtype('double').char,('y'))

    lon = ncfile.createVariable('lon',np.float32,('x'))
    lat = ncfile.createVariable('lat',np.float32,('y'))
    data = ncfile.createVariable(varname,np.float32,('y','x')) # !!!! needed to reverse order for Damon


    # write data to variable.
    lon[:] = lon_in
    lat[:] = lat_in
    data[:] = data_in

    # write global attributes from file
    if len(attrib) > 0:
        af=open(attrib,'r')

        for line in af.readlines():
            attname=line.rsplit('=')[0].replace('\n','').strip()
            attval=line.rsplit('=')[1].replace('\n','').strip()
            
            setattr(ncfile,attname, attval)
            
        af.close()
        


    # close the file.
    ncfile.close()
    print '*** SUCCESS writing '+ncpath


def netcdfappvar(data_in,varname,ncpath):
    ncfile = Dataset(ncpath,'a')

    data = ncfile.createVariable(varname,np.float32,('y','x'))
    data[:] = data_in    
    ncfile.close()

    print '*** SUCCESS appending to '+ncpath
#################################################################




from IPython import embed
from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from rpn.rpn import RPN
from rpn.domains.rotated_lat_lon import RotatedLatLon
import numpy as np

def plotrcm(data,lon,lat,clevs=10,units='',title='',cmap='jet'):
    '''Plot RCMdata on a map where
    - data/lon/lat are 2D arrays of the same size
    - clevs: can be an integer for the number of contours between min and max or a list of values.'''


    # create figure and axes instances
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # create polar stereographic Basemap instance.

    
    custcmap=plt.get_cmap(cmap)        


    if isinstance(clevs, int):
        dmin=np.min(data)
        dmax=np.max(data)
        delbin=(dmax-dmin)/(clevs)
        clevs=np.arange(dmin,dmax+delbin,delbin)
    
    
    if lon[0,0] < -60.:
        domain='QC'
    else:
        domain='EU'

    m = getbasemapfromRPN(domain)

    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    #m.drawstates()
    #m.drawcountries()


    # draw parallels.
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = np.arange(180.,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

    ny = data.shape[0]; nx = data.shape[1]
    x, y = m(lon, lat) # compute map proj coordinates.
    cs = m.contourf(x,y,data,clevs,cmap=custcmap)
    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label('('+units+')')
    # add title
    plt.title(title)
    plt.show()



def getbasemapfromRPN(domain):
    '''Construct basemap object from raw model output file.'''

    vname='MY'
    if domain is 'QC':
        path='/exec/leduc/ClimEx-TINV/kay/pm1950010100_00000000p'
    elif domain is 'EU':
        path='/exec/leduc/ClimEx-TINV/kax/pm1950010100_00000000p'        


    r = RPN(path)
   
    #read a field 
    data = r.get_first_record_for_name(vname)
   
    #get longitudes and latitudes fields
    lons2d, lats2d = r.get_longitudes_and_latitudes_for_the_last_read_rec()

    #get projection parameters
    params = r.get_proj_parameters_for_the_last_read_rec()
   
    #create projection object
    rll = RotatedLatLon(**params)

    #Get the basemap object for the projection and domain defined by the coordinates
    b = rll.get_basemap_object_for_lons_lats(lons2d, lats2d)
        
    return b

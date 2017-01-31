from IPython import embed
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from rpn.rpn import RPN
from rpn.domains.rotated_lat_lon import RotatedLatLon
import numpy as np

def plotrcm(data,lon,lat,clevs=10,units='',title='',cmap='jet',tight=0):
    '''Plot RCMdata on a map where
    - data/lon/lat are 2D arrays of the same size
    - clevs: can be an integer for the number of contours between min and max or a list of values.'''


    # create figure and axes instances
    fig = plt.figure(figsize=(8,8))

    if tight==0:
        ax = fig.add_axes([0.1,0.1,0.8,0.8])

    elif tight==1:
        ax = fig.add_axes([0.,0.06,1,0.9])
    #left=None, bottom=None, right=None, top=None, wspace=None, hspace=None

    # create polar stereographic Basemap instance.

    
    custcmap=plt.get_cmap(cmap)        


    if isinstance(clevs, int):
        dmin=np.min(data)
        dmax=np.max(data)
        delbin=(dmax-dmin)/(clevs)
        clevs=np.arange(dmin,dmax+delbin,delbin)
    

    trim=1
    if lon[0,0] < -60.:
        domain='QC'

        if trim == 0:
            lon0= 83.0
            lon_0= -97.0
            o_lon_p= 180.0
            o_lat_p= 42.5
            llcrnrlon= 260.754
            llcrnrlat= 30.1256
            urcrnrlon= 333.256
            urcrnrlat= 52.8049

        elif trim == 1:
            lon0= 83.0
            lon_0= -97.0
            o_lon_p= 180.0
            o_lat_p= 42.5
            llcrnrlon= 266.204
            llcrnrlat= 34.8294
            urcrnrlon= 322.549
            urcrnrlat= 52.6411

    else:
        domain='EU'

        if trim ==0:
            lon0= -162.0
            lon_0= -342.0
            o_lon_p= 180.0
            o_lat_p= 39.25
            llcrnrlon= 350.983
            llcrnrlat= 24.4208
            urcrnrlon= 51.9427
            urcrnrlat= 66.2037
            
        elif trim ==1:
            lon0= -162.0
            lon_0= -342.0
            o_lon_p= 180.0
            o_lat_p= 39.25
            llcrnrlon= 353.8
            llcrnrlat= 30.428
            urcrnrlon= 37.9136
            urcrnrlat= 63.4463


    m=Basemap(projection="rotpole", lon_0=lon0 - 180, o_lon_p=o_lon_p, o_lat_p=o_lat_p,
              llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
              urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)

    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines(linewidth=2)
    #m.drawstates()
    #m.drawcountries()


    if tight==0:
        # draw parallels.
        parallels = np.arange(0.,90,10.)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=16)
        # draw meridians
        meridians = np.arange(180.,360.,10.)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=16)

    ny = data.shape[0]; nx = data.shape[1]
    x, y = m(lon, lat) # compute map proj coordinates.
    cs = m.contourf(x,y,data,clevs,cmap=custcmap)
    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label('('+units+')')
    # add title
    plt.title(title)

    plt.show()



def getbasemapfromRPN(path,vname='TT'):
    '''Construct basemap object from raw model output file.

    Examples of path:

    /exec/leduc/GRIDS-RPN/kda_dp_280x280
    /exec/leduc/GRIDS-RPN/kba_dp_280x280
    
    '''


    #vname='TT'
    # if domain is 'QC':
    #     #path='/exec/leduc/ClimEx-TINV/kay/pm1950010100_00000000p'
    #     path='/exec/leduc/GRIDS-RPN/kda_dp_280x280'
    # elif domain is 'EU':
    #     #path='/exec/leduc/ClimEx-TINV/kax/pm1950010100_00000000p'
    #     path='/exec/leduc/GRIDS-RPN/kba_dp_280x280'        








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

    # Parameters to call Basemap
    # Basemap(projection="rotpole", lon_0=lon0 - 180, o_lon_p=o_lon_p, o_lat_p=o_lat_p,
    #         llcrnrlon=lons2d[0, 0], llcrnrlat=lats2d[0, 0],
    #         urcrnrlon=lons2d[-1, -1], urcrnrlat=lats2d[-1, -1], **kwargs)

    lon0= rll.get_true_pole_coords_in_rotated_system()[0]
    lon_0=rll.get_basemap_params()['lon_0']
    o_lon_p=rll.get_basemap_params()['o_lon_p']
    o_lat_p=rll.get_basemap_params()['o_lat_p']        
    llcrnrlon=lons2d[0, 0]
    llcrnrlat=lats2d[0, 0]
    urcrnrlon=lons2d[-1, -1]
    urcrnrlat=lats2d[-1, -1]

    print 'lon0=',lon0     
    print 'lon_0=',lon_0    
    print 'o_lon_p=',o_lon_p  
    print 'o_lat_p=',o_lat_p  
    print 'llcrnrlon=',llcrnrlon
    print 'llcrnrlat=',llcrnrlat
    print 'urcrnrlon=',urcrnrlon
    print 'urcrnrlat=',urcrnrlat
    
    return b

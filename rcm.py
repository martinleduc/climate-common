from IPython import embed
from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from rpn.rpn import RPN
from rpn.domains.rotated_lat_lon import RotatedLatLon

def plotrcm(data,lon,lat,clevs,units='',title=''):
    # plot rainfall from NWS using special precipitation
    # colormap used by the NWS, and included in basemap.


    # create figure and axes instances
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # create polar stereographic Basemap instance.

    vv='pr'
    if vv=='tas':
        custcmap=plt.get_cmap('Reds')
    elif vv=='pr':
        custcmap=plt.get_cmap('Blues')
    
    
    if lon[0,0] < -60.:
        domain='QC'



        # My try
        # m = Basemap(projection='lcc', width=5000000, height=5000000,
        #             lat_1=45., lat_2=55, lat_0=48, lon_0=-70.,
        #             resolution = 'l')
        # SEB params
        # m = Basemap(projection='lcc', width=2800000, height=2200000,
        #                       lat_1=45., lat_2=55, lat_0=53, lon_0=-70.,
        #                       resolution = 'l')        

        
    else:
        domain='EU'
        crnrs=[-6.981536865234375,
               29.427963256835938,
               39.913631439208984,
               64.165000915527344]
        m = Basemap(llcrnrlon=crnrs[0],llcrnrlat=crnrs[1],urcrnrlon=crnrs[2],urcrnrlat=crnrs[3],
                    projection='lcc', 
                    lat_1=48.97, lat_2=49.97, lat_0=48.97, lon_0=18.0,
                    resolution = 'l')

    

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
    #lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
    x, y = m(lon, lat) # compute map proj coordinates.
    # draw filled contours.
    #clevs = [0,1,2.5,5,7.5,10,15,20,30]
    #clevs = np.arange(0,36,2)
    #cs = m.contourf(x,y,data,clevs)
    #cs = m.contourf(x,y,data,clevs,cmap=cm.s3pcpn)
    #cs = m.contourf(x,y,data,clevs,cmap=cm.YlOrRd)
    cs = m.contourf(x,y,data,clevs,cmap=custcmap)
    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label('('+units+')')
    # add title
    plt.title(title)
    plt.show()



def getbasemapfromRPN(domain):
    
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

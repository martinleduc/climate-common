#import mpl_toolkits
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys

import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
#from matplotlib.colors import BoundaryNorm
from IPython import embed
################################################################################
# Map plotting utilities
################################################################################


def plotmap(lon,lat,data,units='',lims=[None,None],ax=None,proj='cyl',cmap='jet',cbon=1,crnrs=None,nticks=None):
    '''Plot a map of data (2D) with lon and lat given as arrays. Corners of the map are set according to coordinates or with crnrs = [llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat] and value range by lims[vmin,vmax]
    '''

    f=plt.figure(figsize=(8,4))
    f.subplots_adjust(left=0.02,right=0.88,top=0.95,bottom=0.0)
    if crnrs is None:
        m = Basemap(projection=proj,llcrnrlon=lon[0],llcrnrlat=lat[0],
                    urcrnrlon=lon[-1],urcrnrlat=lat[-1])
    else:
        m = Basemap(projection=proj,llcrnrlon=crnrs[0],llcrnrlat=crnrs[1],
                    urcrnrlon=crnrs[2],urcrnrlat=crnrs[3])
    
    lon_mat, lat_mat = np.meshgrid(lon, lat)
    x, y = m(lon_mat, lat_mat)
    
    m.drawcoastlines(linewidth=0.5)

    vmin=lims[0]
    vmax=lims[1]

    cmap = cmap_discretize(plt.get_cmap(cmap),lims[2])
    #norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cs = m.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap)#, norm=norm)
    #cs = m.pcolormesh(x, y, data, cmap=cmap)

    cb=m.colorbar()
    cb.set_label(units)

    if (lims[0] is not None) and (lims[1] is not None) and (nticks is not None):
        cb.set_ticks(np.linspace(lims[0],lims[1],nticks,endpoint=False))
        

        
    #m.drawmeridians(np.arange(-60, 84, 24), labels = [1,1,1,1]);
    #m.drawparallels(np.arange(-60, 84, 24));

    return cs


######################
# COLORMAP FUNCTIONS #
######################

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    cmap: colormap instance, eg. cm.jet.
    N: number of colors.
    Example
    x = resize(arange(100), (5,100))
    djet = cmap_discretize(cm.jet, 5)
    imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)



    
def custom_div_cmap(numcolors=11, name='custom_div_cmap', cols=['blue', 'white', 'red']):
    """ Create a custom diverging colormap with three colors.
    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()"""

    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list(name=name,
                                             colors = cols,
                                             N=numcolors)
    return cmap


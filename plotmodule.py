from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib
import matplotlib as mpl
import maputil as mu
#import tcrevars
#reload(tcrevars)
###from tcrevars import * !!!!!!!!!!!!!!!!

def dummyplot():
    x=[0,1,2,3,4,5,6]
    y=[6,4,7,1,8,9,7]
    s4 = {'ls':'','marker':'.',"markevery":2,'markersize':20,'markeredgewidth':0}
    #s4 = {'ls':':'}
    
    plt.plot(x,y,**s4)
    plt.show()



    
###### FUNCTIONS related to a Netcdf4 file###########

def plotncdfS(nc,vname,tno):
    ''' Spatial plot of a specific time step (tno) of a given variable (vname) from a netcdf file (nc: tavg type).'''
    # create figure and axes instances
    fig = plt.figure()
    # fig = plt.figure(figsize=(8,8))
    # ax = fig.add_axes([0.1,0.1,0.8,0.8])

    m = Basemap(projection='cyl',lon_0=-90,resolution='c')
    lon2d,lat2d = np.meshgrid(nc.variables['longitude'][:],nc.variables['latitude'][:])


    
    # map feaures
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,90.,45.))
    m.drawmeridians(np.arange(-180.,180.,45.))

    # draw filled contours.
    clevs=myscales(vname)

    if isinstance(clevs,str):
        cs = m.contourf(lon2d,lat2d,nc.variables[vname][tno],latlon='false')
    else:
        cs = m.contourf(lon2d,lat2d,nc.variables[vname][tno],clevs,latlon='false')


    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label(nc.variables[vname].units)

    # add title and show
    plt.title(nc.variables[vname].long_name)

#    ts1=nc.variables['time'][tno]
#    ts2=nc.variables['time'].units

#    print ts1,ts2
#    fig.show()
    return fig


def plotncdfT(nc,vname):
    ''' Plot temporal series of a given variable from a netcdf object (nc: tsi type)'''

    xunits=nc.variables['time'].units
    yunits=nc.variables[vname].units
    
    h=plt.plot(nc.variables['time'][:],nc.variables[vname][:],lw=2)
    #plt.xlim(6000,6099)
    #plt.ylim(13.24,13.31)                 
    plt.ylabel(nc.variables[vname].long_name.split()[-1]+' ('+yunits+')')
    plt.xlabel('time ('+xunits+')')
    plt.title(nc.variables[vname].long_name)


    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    ax.grid(b=True, which='major', linestyle='--')
#    mpl.ticker.ScalarFormatter(useOffset=False)
#    ax.yaxis.set_major_formatter(y_formatter)
#    plt.show()
    return h


    

###### FUNCTIONS related to a Field() object ###########


class Field():
    """ Classe Field:
    attributes: data, sname, lname, fc (conversion factor), units, scale"""

def extTA(ncfile,vname):
    '''Extract time average into a Field object
    ncfile: tavg.nc'''

    nc = Dataset(ncfile)
    ff = Field()
    ff.vname = vname
    ff.data = nc.variables[vname]
    ff.lats = nc.variables['latitude']
    ff.lons = nc.variables['longitude']
    ff.scale = myscales(vname)
    ff.time=nc.variables['time']
    
#    convertfield(ff)
    return ff


def extTS(ncfile,vname):
    '''Extract space averaged time series into a Field object
    ncfile: tsi.nc'''
    nc = Dataset(ncfile)
    ff = Field()
    ff.vname = vname
    ff.data = nc.variables[vname]
    ff.scale = myscales(vname)
    ff.time=nc.variables['time']
    
 #   convertfield(ff)
    return ff

def extTSGP(ncfile,vname):
    '''Extract time series at a specific grid into a Field object
    ncfile: tavg.nc'''
    nc = Dataset(ncfile)
    ff = Field()
    ff.vname = vname
    ff.data = nc.variables[vname]
    ff.scale = myscales(vname)
    ff.data=ff.data[:,80,55]
#    ff.data.units=nc.variables[vname].units
#    ff.data.units='ewrewre'
    ff.time=nc.variables['time']
    
 #   convertfield(ff)
    return ff



# Deprecated - see maputil.py
def plotmap(field,tno):
    ''' Plot a specific time step (tno) of a Field object.'''
    # create figure and axes instances
    fig = plt.figure()
    # fig = plt.figure(figsize=(8,8))
    # ax = fig.add_axes([0.1,0.1,0.8,0.8])

    m = Basemap(projection='cyl',lon_0=-90,resolution='c')
    lon2d,lat2d = np.meshgrid(field.lons,field.lats)

    
    # map feaures
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,90.,45.))
    m.drawmeridians(np.arange(-180.,180.,45.))

    # draw filled contours.
    clevs=field.scale

    if isinstance(clevs,str):
        cs = m.contourf(lon2d,lat2d,field.data[tno],latlon='false')
    else:
        cs = m.contourf(lon2d,lat2d,field.data[tno],clevs,latlon='false')


    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label(field.data.units)

    # add title and show
    plt.title(field.data.long_name)

    #fig.show()
    return fig








######### INFO FUNCTIONS ###########

def myscales(vname):
    if vname is 'A_sat':
        scale=np.arange(-50,40,5)
    elif vname is 'F_precip':
        #scale=[0., 1.,2.,3.,4.,5.]
        scale=np.arange(0,8,1)
    else:
        scale='No scale defined for this variable.'
        
    return scale
    
# def convertfield(ff):
#     if ff.vname is 'F_precip':
#         ff.data=86400*ff.data
#         ff.units='mm/day'




#####################################################
############ Spatial plot from raw data



def plotraw2map(data,tno,scale,title,units,nc,simdir,outdir,fname):
    ''' Plot a specific time step (tno) of a 3-D matrix.'''
    
    # create figure and axes instances
    fig = plt.figure()
    # fig = plt.figure(figsize=(8,8))
    # ax = fig.add_axes([0.1,0.1,0.8,0.8])

    m = Basemap(projection='cyl',lon_0=-90,resolution='c')
    lon2d,lat2d = np.meshgrid(nc.variables['longitude'][:],nc.variables['latitude'][:])

    
    # map feaures
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,90.,45.))
    m.drawmeridians(np.arange(-180.,180.,45.))

    # draw filled contours.
    clevs=np.arange(scale[0],scale[1],scale[2])

#    if isinstance(clevs,str):
#    cs = m.contourf(lon2d,lat2d,data[tno],latlon='false')
#    else:
    cs = m.contourf(lon2d,lat2d,data[tno],clevs,latlon='false')


    # add colorbar.
    cbar = m.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label(units)

    # add title and show
    plt.title(title)

    fig.show()

    if outdir:
        plt.savefig(outdir+fname+'.eps')
        plt.savefig(outdir+fname+'.png',dpi=60)


    return fig





################## Fonction colormap ##############################

# Moved to maputils
# def easycmap(ncol):
    
#     cmap = mpl.cm.jet
#     cstep=cmap.N/float(ncol)
#     cmapdiscrete=[cmap(int(i)) for i in np.arange(0,cmap.N,cstep)]
#     #cmaplist = [cmap(i) for i in range(cmap.N)]
#     #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
#     #bounds = linspace(0,len(my_values),len(my_values)+1)

#     #print np.arange(0,cmap.N-1,cstep)
    
#     return cmapdiscrete


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
        cmap = matplotlib.cm.get_cmap(cmap)
        colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
        colors_rgba = cmap(colors_i)
        indices = np.linspace(0, 1., N+1)
        cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]

    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)



####################### Legend and Color cycle ###########################




def legcolsty(numcol=11,mixid=1):
    '''cx: palette de couleurs
    sx: types de symbole'''

    # c1='b'
    # c2='g'
    # c3='r'
    # c4='m'
    # c5='c'
    # c6='y'
    # c7=(0.6,0.3,0.2)
    # c8=(0.5,0.5,0.5)
    # c9=(1,127./255,0)
    # c10=(0,250./255.,154./255.)
    # c11=(1,20./255,147./255)
    #s4 = {'ls':'--',"marker": "d","markevery":5,'markersize':4,'markeredgewidth':0,'lw':0.}
    #s4 = {'ls':'','marker':'.',"markevery":10,'markersize':10,'markeredgewidth':0}


    c0=(0,0,0)
    s1 = {'ls':'-','lw':3}
    s2 = {'ls':'--','lw':3}
    s3 = {'ls':':','lw':6}

    
    # Mixids
    if mixid==1:
        colv0=mu.easycmap(numcol)
        colv1=colv0
        colv2=[colv0[i] for i in (0,2,4,10)]
        colv=colv0+colv1+colv2+['k']
        lstyv=numcol*[s1]+numcol*[s2]+5*[s3]
    elif mixid==3:
        colv0=mu.easycmap(numcol)
        colv1=colv0
        colv2=[colv0[i] for i in (0,2,4,10)]
        colv=colv0+colv1+colv2+['k']        
        lstyv=numcol*[s1]+numcol*[s2]+5*[s1]        
    elif mixid==2:
        offsize=8
        colv0=mu.easycmap(numcol)
        colv1=colv0
        colv2=mu.easycmap(offsize)              #8 param!
        colv=colv0+colv1+colv2
        lstyv=numcol*[s1]+numcol*[s1]+offsize*[s2]
    elif mixid==4:
        # RCP only, numcol = num PE and CE
        colv0=mu.easycmap(numcol)
        colv1=colv0
        colv2=[colv0[i] for i in (0,2,4,10)]
        colv=colv2+colv2
        lstyv=4*[s1]+4*[s2]        

    
    
    # ver pour 11+4+1+4
    # colv0=easycmap(11)
    # colv1=[colv0[i] for i in (0,2,4,10)]+['k']+[colv0[i] for i in (0,2,4,10)]
    # colv=colv0+colv1
    # lstyv=11*[s1]+5*[s2]+4*[s3]

    # pour version 1 du papier
    # colv0=mu.easycmap(numcol)
    # colv1=mu.easycmap(numcol)
    # colv2=[colv0[i] for i in (0,2,4,10)]
    # colv=colv0+colv1+colv2+[c0]+[c0]
    # lstyv=numcol*[s1]+numcol*[s2]+4*[s4]+2*[s4]



    #+4*[s4]+2*[s4]
    

    
    # colv1=[colv0[i] for i in (0,2,4,10)]+['k']+[colv0[i] for i in (0,2,4,10)]

    
    # colv=colv0+colv1
    

    
    # colv1=[colv0[i] for i in (0,2,4,10)]+['k']+[colv0[i] for i in (0,2,4,10)]
    # colv=colv0+colv1
    # lstyv=11*[s1]+5*[s2]+4*[s3]


    #lstyv=lstyv[mixid:]

    # General Graphics parameters
    #mpl.rcParams['axes.color_cycle']=colv[mixid:]
    mpl.rcParams['axes.color_cycle']=colv
    mpl.rcParams['legend.frameon']='False'
    mpl.rcParams.update({'font.size': 24})
    #mpl.rcParams.update({'legend.fontsize': 11})
    mpl.rcParams.update({'legend.fontsize': 18})
    mpl.rcParams.update({'grid.linewidth': 0.5})

    
    return lstyv



# def legendcolors(bid=-1):


#     c1='b'
#     c2='g'
#     c3='r'
#     c4='c'
#     c5='m'
#     c6='y'
#     c7=(0.6,0.3,0.2)
#     c8=(0.5,0.5,0.5)
#     c9=(1,127./255,0)
#     c10=(0,250./255.,154./255.)
#     c11=(1,20./255,147./255)
    
    
#     #mpl.rcParams['axes.color_cycle']=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11]
#     mpl.rcParams['axes.color_cycle']=[c1,c2,c3,c4,c5,c6,c2,c4,c6,c7,c11]
    
#     print str(mpl.rcParams['axes.color_cycle']) #mpl.rcParams['axes.color_cycle']=['b','g','r','c','m','y',(0.6,0.3,0.2),(0.5,0.5,0.5),(1,127./255,0),(0,250./255.,	154./255.)]

#     if bid >= 0:
#         #mpl.rcParams['axes.color_cycle'].pop(bid)
#         mpl.rcParams['axes.color_cycle'].insert(bid,'k')
#     else:
#         mpl.rcParams['axes.color_cycle'].append('k')
#     print str(mpl.rcParams['axes.color_cycle'])

#     mpl.rcParams['legend.frameon']='False'
            
        
#     #mpl.rcParams['axes.color_cycle'].insert(bid,mpl.rcParams['axes.color_cycle'].pop(-1))


def testcolorcycle(nc):

    #legendcolors()
    legcolsty()
    nx=100
    x=np.arange(0,nx)
    plt.figure()

    tlist=[]
    
    for i in np.arange(0,nc):
        #y=i*np.random.randn(nx)
        y=(i+0.1)*x
        plt.plot(x,y,lw=2)
        tlist.append(str(i))

    plt.legend(tlist,2)
    plt.show()







################## Fonction diagnostiques ########################


def diagemi(ncfile):
    
    nc = Dataset(ncfile)

    nt=nc.variables['time'].size
    diagem=np.zeros(nt)
    
    for tt in np.arange(1,nt):
        dA=nc.variables['A_totcarb'][tt] - nc.variables['A_totcarb'][tt-1] 
        dL=nc.variables['L_totcarb'][tt] - nc.variables['L_totcarb'][tt-1]
        dO=nc.variables['O_totcarb'][tt] - nc.variables['O_totcarb'][tt-1]
        diagem[tt]=dA+dL+dO


    return diagem
        
    



#if __name__ == "__main__":
#    import sys
#    plottime(int(sys.argv[1])).show()
#    raw_input('Press enter.')

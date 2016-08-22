from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#import cmapdemo as cmapd



# Need to test with uvic data
# from netCDF4 import Dataset
# b=Dataset("/home/leduc/data/CCR/C1.6/tavg.01851.01.01.nc")


def plotmesh(lons,lats,data,tno=None,var='',sea='',diag=0,units='',title='',pcpch='',ax=None,cbon=1,fs=40):
    if tno is None:
        zdata=data.T
    else:
        zdata=data[tno][0].T



    nx = lons.shape[0]
    ny = lats.shape[0]

    dx = (lons[-1] - lons[0]) / float(nx - 1)
    dy = (lats[-1] - lats[0]) / float(ny - 1)

    #data = np.random.randn(nx, ny)

    lons_p1 = list(lons) + [lons[-1] + dx, ]
    lons_m1 = [lons[0] - dx] + list(lons)
    lons_corn = [(a + b) / 2.0 for a, b in zip(lons_p1, lons_m1) ]

    lats_p1 = list(lats) + [lats[-1] + dy, ]
    lats_m1 = [lats[0] - dy] + list(lats)
    lats_corn = [(a + b) / 2.0 for a, b in zip(lats_p1, lats_m1) ]

    lat2d, lon2d = np.meshgrid(lats_corn, lons_corn)


    clevs,units = mycolorspaces(var,sea,diag=diag)
    
    
    b = Basemap(llcrnrlon=-145,llcrnrlat=20,urcrnrlon=-50,urcrnrlat=72)
    #b = Basemap()
    x, y  = b(lon2d, lat2d)



    if clevs.size == 1:
        cmap = cm.get_cmap("jet", 10)
        cs=b.pcolormesh(x, y, zdata)
        if cbon: cb=plt.colorbar(cs,shrink=0.8)
    else:
        cmap = cm.get_cmap("jet", clevs.size-1)
        cs=b.pcolormesh(x, y, zdata,vmin=clevs[0],vmax=clevs[-1],cmap=cmap)
        if cbon: cb=plt.colorbar(cs,ticks=clevs,shrink=1.)
        
    if cbon: cb.set_label(units)

    plt.title(title,fontsize=fs)

    b.drawcoastlines()
    #b.drawmeridians(np.arange(-60, 84, 24), labels = [1,1,1,1]);
    #b.drawparallels(np.arange(-60, 84, 24));

    plt.tight_layout()    
   

    return cs,clevs


#########################################
def plotmap(lo,la,data,tno,var='',sea='',diag=0,units='',title='',pcpch=''):
    ''' Plot a specific time step (tno) of data(4D) with lat and lon
    --- Should be modified for data(2D)
    
    The shape of lat and lon should be of the form (L,)
    '''
    #fig = plt.figure()
    #fig.canvas.manager.window.Move(100,400)

    #ax = fig.add_axes([0.1,0.1,0.8,0.8])

    # Map
    #m=Basemap(llcrnrlon=217.5,llcrnrlat=22,urcrnrlon=307.5,urcrnrlat=70,projection='mill')



    
    # If data is 4D
    zdata = data[tno][0]
    
    lo, la = np.meshgrid(lo, la)


    
    clevs,units = mycolorspaces(var,sea,diag=diag)

    print '+++++++++'
    print var
    print sea
    print diag
    print clevs
    
    m = Basemap(llcrnrlon=-145,llcrnrlat=20,urcrnrlon=-50,urcrnrlat=72)
    X,Y = m(lo, la)
    
    cs = m.contourf(X, Y, zdata,clevs)
    cb=m.colorbar(cs)
    cb.set_label(units)
    #fig.colorbar(cs, ax=axs, format="%.2f")

    plt.title(title)


    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,90.,45.))
    m.drawmeridians(np.arange(-180.,180.,45.))

    
    #plt.ion()
    #plt.show()

    #return fig




def easycmap(ncol):
    
    cmap = mpl.cm.jet
    cstep=cmap.N/float(ncol)
    cmapdiscrete=[cmap(int(i)) for i in np.arange(0,cmap.N,cstep)]
    #cmaplist = [cmap(i) for i in range(cmap.N)]
    #cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    #bounds = linspace(0,len(my_values),len(my_values)+1)

    #print np.arange(0,cmap.N-1,cstep)
    
    return cmapdiscrete



def savecolorbar(cs,clevs,units,output=''):
    # cs is a mapable object, clevs the levels

    
    plt.ion()


    a = np.array([[0,1]])
    plt.figure(figsize=(3, 20))
    img = plt.imshow(a, cmap="jet")
    plt.gca().set_visible(False)
    
    cax = plt.axes([0.05, 0.05, 0.2, 0.9])
    cb=plt.colorbar(cs,cax=cax,ticks=clevs)
    cb.set_label(units,fontsize=40)
    cb.ax.tick_params(labelsize=40)
    if output:
        plt.savefig(output)



    
    #fcb=plt.figure(figsize=(6,10))
    # fcb=plt.figure(figsize=(6,10))
    # ax=plt.subplot(111)
    # #plt.gca().set_visible(False)
    # plt.colorbar(cs,cax=ax,ticks=clevs)

    # #ax.colorbar(cs,ticks=clevs,pad=0.5)
    # # cb=plt.colorbar(cs,ticks=clevs,pad=0.5)
    # cb.set_label(units)


    
    # plt.gca().set_visible(False)
    # cb=plt.colorbar(cs,ticks=clevs,pad=0.5)
    # cb.set_label(units)

    
    # Colorbar
    # import matplotlib as mpl
    # import matplotlib.pyplot as plt
    # import numpy as np

    #cmapd.discrete_cmap(11)

    #a=np.array([[0,1,2,3,4,5,6]])
    # a=np.array([clevs])
    # fcb=plt.figure(figsize=(3.4,9))
    # #img = plt.imshow(a, cmap="indexed")
    # img = plt.imshow(a)
    # plt.gca().set_visible(False)

    # cax = plt.axes([0.4, 0.05, 0.1, 0.9])
    # cb=plt.colorbar(cax=cax)
    # cb.set_ticks([])
    # cb.set_ticklabels([])


    plt.show()

    #fcb.savefig(figdir+"colorbar.eps")







# MUST DISAPEAR, should be project dependent, give as an argument to mesh
def mycolorspaces(var,sea,diag=0,pcpch=''):

    #Default values
    clevs=np.array([10])
    units='(no units)'

    if var is 'tas':
        units='($^\circ$C)'   
    elif var is 'pr':        
        units='(mm/day)'
        if pcpch is 'perc':
            units='(%)'


    if diag==0: # no delta no diff
        if var=='tas':
            clevs=np.linspace(-5, 35, 10)
        elif var=='tas' and sea=='h':
            clevs=np.linspace(-5, 35, 10)

        elif var=='pr' and sea=='e':
            clevs=np.linspace(0, 11, 12)
            
    elif diag is 'del':
        if var=='tas' and (sea=='e' or sea=='JJA'):        
            #clevs=np.linspace(-2, 9, 12)
            clevs=np.linspace(-1, 8, 13)
        elif var=='tas' and sea=='h':
            clevs=np.linspace(-1, 13, 15)
        elif var=='tas' and sea=='y':
            clevs=np.linspace(0, 6, 11)            

        elif var=='pr' and sea=='e':
            clevs=np.linspace(-1, 1, 11)
        elif var=='pr' and sea=='h':
            clevs=np.linspace(-1, 1.5, 11)
        elif var=='pr' and sea=='y':
            clevs=np.linspace(-0.5, 0.5, 11)

                
    elif diag is 'deldiff':
        if var=='tas' and (sea=='e' or sea=='JJA'):
            clevs=np.linspace(-4, 4, 9)
        elif var=='tas' and sea=='h':
            clevs=np.linspace(-10, 10, 11)
        if var=='pr' and sea=='e':
            clevs=np.linspace(-1.5, 1.5, 13)
        if var=='pr' and sea=='h':
            clevs=np.linspace(-1.5, 1.5, 13)            
            
    elif diag is 'iv':
        if var=='tas' and (sea=='e' or sea=='JJA'):
            clevs=np.linspace(0, 1, 9)
        elif var=='tas' and sea=='h':
            clevs=np.linspace(0, 3, 13)
            
        if var=='pr' and sea=='e':
            clevs=np.linspace(0, 0.5, 11)
        if var=='pr' and sea=='h':
            clevs=np.linspace(0, 1, 11)                        
    elif diag is 'delrdiff':
        if var=='tas' and sea=='e':
            clevs=np.linspace(-100,100, 11)
            units='(%)'

    elif diag is 'norm':
            clevs=np.linspace(0,1, 11)
            units=''

    elif diag is 'uhalfnorm':
            clevs=np.linspace(50.,100., 11)
            #clevs=np.linspace(0,0.5, 11)
            units='($100\\times N\'/N$)'
    elif diag is 'sig2090':
        if var=='tas' and sea=='y':            
            clevs=np.linspace(0,2, 11)
            #clevs=np.linspace(0,0.5, 11)

        elif var=='pr' and sea=='y':            
            clevs=np.linspace(0,0.5, 11)
            #clevs=np.linspace(0,0.5, 11)

    elif diag is 'sig2noise':
        if var=='tas' and sea=='h':            
            clevs=[0,20]
        elif var=='tas' and sea=='e':            
            clevs=[0,20]
        elif var=='tas' and sea=='y':            
            clevs=[0,20]                        
        elif var=='pr' and sea=='h':            
            clevs=[-10,10]
        elif var=='pr' and sea=='e':            
            clevs=[-10,10]

    return clevs,units


    

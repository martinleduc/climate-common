from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
import rcm as rcm
reload(rcm)

sourcedir='/exec/leduc/ClimEx/'
run='kda'
vv='tas'
tstamp='195507'

file=sourcedir+run+'/series/'+tstamp+'/'+vv+'_'+run+'_'+tstamp+'_se.nc'


ff=Dataset(file)

data=ff.variables[vv][0,:,:]-273.15
units='$^\circ$C'

lat=ff.variables['lat'][:]
lon=ff.variables['lon'][:]

title=run
#clevs=np.arange(0,50,10)
rcm.plotrcm(data,lon,lat,clevs=20,units=units,title=title)



# Plot one grit point

yy,xx=100,100
data2=data*0
data2[yy,xx]=1
rcm.plotrcm(data2,lon,lat,title=title)



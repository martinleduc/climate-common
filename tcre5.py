from netCDF4 import Dataset
import numpy as np


#import classModel.Ensemble as Ens
import ensemble.Model as Ens
import ensemble.gridutils as gu


def getlsmask():
    '''Return numeric land-sea mask with land and sea areas represented by 0 (mask off) and 1 (mask on) respectively'''
    
    nc=Dataset('/dmf2/scenario/external_data/CMIP5/CCCMA/CanESM2/1pctCO2/fx/atmos/r0i0p0/areacello/areacello_fx_CanESM2_1pctCO2_r0i0p0.nc')
    
    data=nc.variables['areacello'][:].data
    lon=nc.variables['lon'][:]
    lat=nc.variables['lat'][:]

    lon_dest,lat_dest,data=gu.trlonlatdatacmip5(lon,lat,data)
    #data, lon_dest, lat_dest=Ens.interpolateto(data,lon_dest,lat_dest,[1,1])
    data, lon_dest, lat_dest=gu.interpolateto(data,lon_dest,lat_dest,'CanESM')    

    mask=np.round(data/np.max(data))

    mask=1-mask
    return mask




models=['BNU-ESM',
        'CanESM2',
        'CESM1-BGC',
        'HadGEM2-ES',
        'inmcm4',
        'IPSL-CM5A-LR',
        'IPSL-CM5A-MR',
        'IPSL-CM5B-LR',
        'MIROC-ESM',
        'MPI-ESM-LR',
        'MPI-ESM-MR',
        'NorESM1-ME'
        ]

giorgiregions=['AUS',
               'AMZ',
               'SSA',
               'CAM',
               'WNA',
               'CNA',
               'ENA',
               'ALA',
               'GRL',
               'MED',
               'NEU',
               'WAF',
               'EAF',
               'SAF',
               'SAH',
               'SEA',
               'EAS',
               'SAS',
               'CAS',
               'TIB',
               'NAS',
               ]

mylandregions=['Qc','ANT']

oceanregions=['SPO',
              'NPO',
              'SAO',
              'NAO',
              'IO',
              'AO',
              'SO']


arearegions=['ARC']



lsmask=getlsmask()



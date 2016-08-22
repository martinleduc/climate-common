# Module for manipulating netcdf files
# NOTE: The UVic model nc files are of the type NETCDF3_CLASSIC
from netCDF4 import Dataset
import os
import numpy as np
import plotmodule as pm
import sys as sys
#import netCDF4
#from netCDF4 import Dataset as NCDFF



def cpcdf(infile,outfile):
    # Copy input file to output
    os.system("cp "+infile+" "+outfile)

    


def genco2ccn(infile,outfile,exp):
    ''' Generate 1%-increase of CO2 ccn and make nc file to force the model.'''

    
    cpcdf(infile,outfile)


    if exp=='4x':
        # y1=850.5                                # starting year for 1% CO2 increases
        # y2=y1+140.                               # starting year for stabilisation
        # y3=2005.5                                # end of forcing data
        y1=1850.5                                # starting year for 1% CO2 increases
        y2=y1+140.                               # starting year for stabilisation
        y3=2400.5                                # end of forcing data
    elif exp=='2x':
        # y1=850.5                                # starting year for 1% CO2 increases
        # y2=y1+70                               # starting year for stabilisation
        # y3=2005.5                                # end of forcing data
        y1=1850.5                                # starting year for 1% CO2 increases
        y2=y1+70                               # starting year for stabilisation
        y3=2400.5                                # end of forcing data




    # Open and modify the output

    nc = Dataset(outfile, 'a')
    years=nc.variables['time'][:]

    y1id=np.where(years==y1)
    y2id=np.where(years==y2)
    y3id=np.where(years==y3)


    print years[y1id]
    print years[y2id[0][0]]
    print years[y3id[0][0]]

    
    for tt in range(y1id[0][0]+1,y2id[0][0]):
        nc.variables['A_co2'][tt] = 1.01*nc.variables['A_co2'][tt-1]

    for tt in range(y2id[0][0],y3id[0][0]+1):
        nc.variables['A_co2'][tt] = nc.variables['A_co2'][tt-1]
        
    nc.close()



def genpulsemit(outdir,timeinfo,pulseTt):
    ''' Generate pulse of CO2 emissions and make nc file to force the model.

    Typical values:
    
    pulseTt=1.
    t1=1850.
    t2=1860.
    t3=2350.
    '''
    # PARAMETERS
    emit0=17656.012               # reference level of CO2 emissions (1850) (in kg/s = 0.557 Pg/yr)
    emit1=0.
    pmagn= pulseTt*31709791.9838              # in Kg/s -> N x 1Tt/yr
    var='F_co2emit'

    # Extraire informations temporelles
    t1=timeinfo[0]
    t2=timeinfo[1]
    t3=timeinfo[2]

    outfile=var+'-pulse'+str(int(round(t2)))+'-'+str(pulseTt)+'Tt'+'.nc'
    print outdir+outfile

    # Calculer axe temporel
    timeaxis=np.arange(t1,t3+1,1)
    nt=timeaxis.size

    pulsetimeid=np.where(timeaxis==t2)[0][0]

    
    # Create NETCDF file
    os.path.exists(outdir+outfile) and os.remove(outdir+outfile)
    nc = Dataset(outdir+outfile, 'w',format='NETCDF3_CLASSIC')

    nc.createDimension('time',nt)
    time=nc.createVariable('time',np.dtype('float64').char,('time'))
    time[:]=timeaxis[:]
    time.axis='T'
    time.long_name='time'
    time.standard_name='time'
    time.units='year'
    #time.calendar="noleap"

    data1=nc.createVariable('F_co2efuel',np.dtype('float64').char,('time'))
    data1[:]=np.zeros(nt)+emit0
    data1[-(nt-pulsetimeid):]=np.zeros(nt-pulsetimeid)+emit1
    data1[pulsetimeid]=pmagn
    data1.long_name='co2 emissions from fossil fuels'
    data1.units='kg s-1'

    
    data2=nc.createVariable('F_co2eland',np.dtype('float64').char,('time'))
    data2[:]=np.zeros(nt)
    data2.long_name='co2 emissions from land changes'
    data2.units='kg s-1'

    
    data3=nc.createVariable('F_co2emit',np.dtype('float64').char,('time'))
    data3[:]=data1[:]+data2[:]
    data3.long_name='co2 emissions'
    data3.units='kg s-1'

    nc.close()



# def genfunemit(outdir,timeinfo,cumTt,typef):
#     ''' Generate linear function for determining CO2 cumulative emissions and make nc file to force the model with.

#     Typical values:
    
#     '''

#     # PARAMETERS
#     emit0=17656.012               # reference level of CO2 emissions (1850) (in kg/s = 0.557 Pg/yr)
#     cum=cumTt*31709791.9838              # in Kg/s -> N x 1Tt/yr
#     var='F_co2emit'


#     # Extraire informations temporelles
#     t1=timeinfo[0]
#     t2=timeinfo[1]
#     t3=timeinfo[2]

#     outfile=var+"-"+typef+str(int(round(t2)))+'-'+str(cumTt)+'Tt'+'.nc'
#     #outfile=var+'.nc'

#     # Calculer axe temporel
#     timeaxis=np.arange(t1,t3+1,1)
#     nt=timeaxis.size

#     t2timeid=np.where(timeaxis==t2)[0][0]

    
#     # Create NETCDF file
#     os.path.exists(outdir+outfile) and os.remove(outdir+outfile)
#     nc = Dataset(outdir+outfile, 'w',format='NETCDF3_CLASSIC')

#     nc.createDimension('time',nt)
#     time=nc.createVariable('time',np.dtype('float64').char,('time'))
#     time[:]=timeaxis[:]
#     time.axis='T'
#     time.long_name='time'
#     time.standard_name='time'
#     time.units='year'
#     #time.calendar="noleap"

#     data1=nc.createVariable('F_co2efuel',np.dtype('float64').char,('time'))
#     data1.long_name='co2 emissions from fossil fuels'
#     data1.units='kg s-1'

#     if typef=="linear":
#         emit1=cum/(nt-t2timeid)
#         data1[:]=np.zeros(nt)+emit0
#         data1[-(nt-t2timeid):]=np.zeros(nt-t2timeid)+emit1
        
    
#     data2=nc.createVariable('F_co2eland',np.dtype('float64').char,('time'))
#     data2[:]=np.zeros(nt)
#     data2.long_name='co2 emissions from land changes'
#     data2.units='kg s-1'

    
#     data3=nc.createVariable('F_co2emit',np.dtype('float64').char,('time'))
#     data3[:]=data1[:]+data2[:]
#     data3.long_name='co2 emissions'
#     data3.units='kg s-1'

#     nc.close()

def genfunemit(outdir,timeinfo,cumTt,typef):
    ''' Generate linear constant rate emission scenarios (t1,t2,t3) with the ability to switchoff (t1,t2,toff,t3) the emissions at some point. Generate an nc file to force the model with.


    '''

    # PARAMETERS
    emit0=17656.012               # reference level of CO2 emissions (1850) (in kg/s = 0.557 Pg/yr)
    cum=cumTt*31709791.9838              # in Kg/s -> N x 1Tt/yr
    var='F_co2emit'


    # Extraire informations temporelles
    if typef=="linear":
        t1=timeinfo[0]
        t2=timeinfo[1]
        t3=timeinfo[2]
        outfile=var+"-"+typef+str(int(round(t2)))+'-'+str(cumTt)+'Tt'+'.nc'        
    elif typef=="sawteeth":
        t1=timeinfo[0]
        t2=timeinfo[1]
        toff=timeinfo[2]
        t3=timeinfo[3]
        outfile=var+"-"+typef+str(int(round(t2)))+'-'+str(int(round(toff)))+'-'+str(cumTt)+'Tt'+'.nc'

        

    #outfile=var+'.nc'

    # Calculer axe temporel
    timeaxis=np.arange(t1,t3+1,1)
    nt=timeaxis.size

    t2timeid=np.where(timeaxis==t2)[0][0]
    tofftimeid=np.where(timeaxis==toff)[0][0]    
    
    # Create NETCDF file
    os.path.exists(outdir+outfile) and os.remove(outdir+outfile)
    nc = Dataset(outdir+outfile, 'w',format='NETCDF3_CLASSIC')

    nc.createDimension('time',nt)
    time=nc.createVariable('time',np.dtype('float64').char,('time'))
    time[:]=timeaxis[:]
    time.axis='T'
    time.long_name='time'
    time.standard_name='time'
    time.units='year'
    #time.calendar="noleap"

    data1=nc.createVariable('F_co2efuel',np.dtype('float64').char,('time'))
    data1.long_name='co2 emissions from fossil fuels'
    data1.units='kg s-1'

    if typef=="linear":
        emit1=cum/(nt-t2timeid)
        data1[:]=np.zeros(nt)+emit0
        #data1[-(nt-t2timeid):]=np.zeros(nt-t2timeid)+emit1
        data1[t2timeid:]=np.zeros(nt-t2timeid)+emit1
    if typef=="sawteeth":
        emit1=cum/(nt-t2timeid)
        data1[:]=np.zeros(nt)+emit0
        data1[t2timeid:tofftimeid]=np.zeros(tofftimeid-t2timeid)+emit1
        data1[tofftimeid:]=0
        
    
    data2=nc.createVariable('F_co2eland',np.dtype('float64').char,('time'))
    data2[:]=np.zeros(nt)
    data2.long_name='co2 emissions from land changes'
    data2.units='kg s-1'

    
    data3=nc.createVariable('F_co2emit',np.dtype('float64').char,('time'))
    data3[:]=data1[:]+data2[:]
    data3.long_name='co2 emissions'
    data3.units='kg s-1'

    nc.close()



def gen1pctccn(outdir,timeinfo,ccn0,timeco2):
    ''' Generate mathematical function for determining CO2 concentrations and make nc file to force the model with.

    A_1pct2xCO2-1850.nc'
    
    '''

    # PARAMETERS
    #emit0=17656.012               # reference level of CO2 emissions (1850) (in kg/s = 0.557 Pg/yr)
    #cum=cumTt*31709791.9838              # in Kg/s -> N x 1Tt/yr
    var='A_co2'


    # Extraire informations temporelles
    t1=timeinfo[0]
    t2=timeinfo[1]
    t3=timeinfo[2]

    outfile=var+"-1pct-"+str(timeco2)+"xCO2"+'-1850.nc'

    # Calculer axe temporel
    timeaxis=np.arange(t1,t3+1,1)
    nt=timeaxis.size

    t2timeid=np.where(timeaxis==t2)[0][0]

    
    # Create NETCDF file
    os.path.exists(outdir+outfile) and os.remove(outdir+outfile)
    nc = Dataset(outdir+outfile, 'w',format='NETCDF3_CLASSIC')

    nc.createDimension('time',nt)
    time=nc.createVariable('time',np.dtype('float64').char,('time'))
    time[:]=timeaxis[:]
    time.axis='T'
    time.long_name='time'
    time.standard_name='time'
    time.units='year'
    #time.calendar="noleap"
    
    data1=nc.createVariable('A_co2',np.dtype('float64').char,('time'))
    data1.long_name='CO2 concentration'
    data1.units='ppmv'

    data1[:]=np.zeros(nt)+ccn0
    for tt in np.arange(t2timeid,nt):
        ccn1=1.01*data1[tt-1]

        if ccn1 >= timeco2*data1[0]:
            data1[tt]=data1[tt-1]
        else:
            data1[tt]=ccn1

    nc.close()

    
   

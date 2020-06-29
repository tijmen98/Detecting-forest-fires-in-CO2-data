# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:23:08 2020

@author: tijme
"""
#import modules
from netCDF4 import Dataset
import numpy as np

#%%

#data directories
loc= 'E:\Data'
readname='\data4.nc'
writename='\data4new.nc'

#%%

#import data
Data= Dataset(loc+readname,'r')
lat= Data.variables['latitude'][:]
lon= Data.variables['longitude'][:]
time=Data.variables['time'][:]

print('data imported')

#create dataset
newdata= Dataset(loc+writename,'w',format="NETCDF4")
newdata= Dataset(loc+writename,'a')
#create dimensions
newlat=newdata.createDimension('latitude',len(lat))
newlon=newdata.createDimension('longitude',len(lon))
newtime=newdata.createDimension('time',None)
#create variables
latitudes = newdata.createVariable("latitude","f4",('latitude',))
longitudes = newdata.createVariable("longitude","f4",('longitude',))
times= newdata.createVariable("time","f8",('time',))
Windspeeds = newdata.createVariable("windspeed","f8",("time","latitude","longitude",))
Directions = newdata.createVariable("winddirection","f8",("time","latitude","longitude",))
#assign data to dimensions
latitudes[:]=lat[:]
longitudes[:]=lon[:]
times[:]=time[:]
#close new dataset
newdata.close()

print('Netcdf made')

#%%
#Def programs

#direction
def direction(V,U):
    A = 0
    if V>0:
        if U>0:
            #North-east quadrant
            A =((np.arctan(abs(U)/abs(V))*180)/np.pi)
        else:
            #North-west quadrant
            A = 360 - ((np.arctan(abs(U)/abs(V))*180)/np.pi)
    else:
        if U>0:
            #South-east quadrant
            A = 180 - ((np.arctan(abs(U)/abs(V))*180)/np.pi)
        else:
            #South-west quadrant
            A = 180 + ((np.arctan(abs(U)/abs(V))*180)/np.pi)
    return(A)
                        
                        
                
            

#windspeed
def speed(V,U):
    Windspeed = np.sqrt(V**2+U**2)
    if Windspeed > 1000:
        Windspeed = -1
    
    return(Windspeed)

#%%
    

#%%

lentime=len(time)
lenlat=len(lat)
lenlon=len(lon)

V=0
U=0

singlespeed = np.zeros([lenlat,lenlon])
singledirection = np.zeros([lenlat,lenlon])
for t in np.arange(0,(lentime),1):
    v=Data.variables['v10'][t,:,:] #north vector
    V=np.array(v)
    u=Data.variables['u10'][t,:,:] #north vector
    U=np.array(u)
    print(t, 'from',lentime)
    for x in np.arange(0,(lenlat),1):
        for y in np.arange(0,(lenlon),1):
            singlespeed[x,y]=speed(V[x,y],U[x,y])
            singledirection[x,y]=direction(V[x,y],U[x,y])
    newdata= Dataset(loc+writename,'a')    
    Windspeeds=newdata.variables['windspeed']
    Directions=newdata.variables['winddirection']
    Windspeeds[t,:,:]=singlespeed
    Directions[t,:,:]=singledirection
    newdata.close()
    
            
            
            
            
            
print('Windspeed and direction calculated')        



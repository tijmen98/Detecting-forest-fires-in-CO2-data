# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:23:08 2020

@author: tijme
"""
#import modules
from netCDF4 import Dataset
import numpy as np
import math
import shapefile

#%%

#data directories read
readloc= 'E:\Data\Winddata'
readname=['\data1new.nc','\data2new.nc','\data3new.nc','\data4new.nc']

#data directories write

writeloc='E:\Data\shapefiles_wind'
writename = '\windpaden25-8-1st'



#create shapefile

file=shapefile.Writer(writeloc+writename)
file.field('afwijking','F',decimal=20)
#file.field('time','F',decimal=20)
 

#%%

"import data"

#import timelenghts
datalens=np.zeros(4)
for step in range(4):
    Data= Dataset(readloc+readname[step],'r')
    time=Data.variables['time'][:]
    datalens[step]=len(time)
    Data.close()

print('data imported')

#%%
"define functions"
    
def newlat(windspeed,bearing,lat,lon):
    R = 6378.1 #Radius of the Earth
    brng = math.radians(bearing) #Bearing is 90 degrees converted to radians.
    d = windspeed*3.6 #Distance in km
    
    lat1 = math.radians(lat) #Current lat point converted to radians
    lon1 = math.radians(lon) #Current long point converted to radians
    
    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
         math.cos(lat1)*math.sin(d/R)*math.cos(brng))
    
    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
                 math.cos(d/R)-math.sin(lat1)*math.sin(lat2))
    
    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)
                
    return(lat2)

def newlon(windspeed,bearing,lat,lon):
    R = 6378.1 #Radius of the Earth
    brng = math.radians(bearing) #Bearing is 90 degrees converted to radians.
    d = windspeed*3.6 #Distance in km
    
    lat1 = math.radians(lat) #Current lat point converted to radians
    lon1 = math.radians(lon) #Current long point converted to radians
    
    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
         math.cos(lat1)*math.sin(d/R)*math.cos(brng))
    
    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
                 math.cos(d/R)-math.sin(lat1)*math.sin(lat2))
    
    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)
            
    return(lon2)    

print('functions defined')

#%%

"initial conditions"
#time in hours
t_start=588
dt=1
t_end=t_start+24
#starting coordinates
startlatitude=71.19    #degrees north
startlongitude=144.72       #degrees east

#afwijking
afwijking = 0.1
#%%


afwijkingarray = [1]
afwijkingarray[1:] = np.random.normal(1,afwijking,100)

for deviate in afwijkingarray: 
    latitude=[startlatitude]
    longitude=[startlongitude]
    timelist=[]
    t=t_start
    while t<t_end:
        if t < datalens[0]:
            Data= Dataset(readloc+readname[0],'r')
            new_t=t
        else:
            if t < datalens[1]+datalens[0]:
                Data= Dataset(readloc+readname[1],'r')
                new_t=t-datalens[0]
            else:
                if t < datalens[2]+datalens[1]+datalens[0]:
                    Data= Dataset(readloc+readname[2],'r')
                    new_t=t-datalens[0]-datalens[1]
                else:
                    if t < datalens[3]+datalens[2]+datalens[1]+datalens[0]:
                        Data= Dataset(readloc+readname[3],'r')
                        new_t=t-datalens[2]-datalens[1]-datalens[0]
        
        lat=latitude[-1]
        lon=longitude[-1]
        windspeed=Data.variables['windspeed'][new_t,int((90-lat)*10),int(lon*10)]
        bearing=Data.variables['winddirection'][new_t,int((90-lat)*10),int(lon*10)]
        if windspeed>0:
            latitude.append(newlat((windspeed*deviate),(bearing*deviate),lat,lon))
            longitude.append(newlon((windspeed*deviate),(bearing*deviate),lat,lon))
        else:
            latitude.append(lat)
            longitude.append(lon)
        timelist.append(t)
        Data.close()
        t=t+1
    vector=[]
    for step in range(t_end-t_start):
        vector.append([longitude[step],latitude[step],timelist[step]])
        #file.point(longitude[step],latitude[step])
        #file.record(deviate,timelist[step])
    file.line([vector])
    file.record(deviate)
    
file.close()
    
    
        
#%%
    
"write data"
""""
file=shapefile.Writer(loc+'\windpaden')
file.field('1e run','C')
file.line(allvectors)
file.record('lijn')
file.close()
 """   

    
    
    
    


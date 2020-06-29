# -*- coding: utf-8 -*-
"""
Created on Thu May 21 15:19:13 2020

@author: Tijmen van der Drift
"""


import numpy as np
import math
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#%%
"data directories"

#data directories read
readloc= 'E:\Data\Winddata'
readname=['\data1new.nc','\data2new.nc','\data3new.nc','\data4new.nc']

#import timelenghts
datalens=np.zeros(4)
for step in range(4):
    Data= Dataset(readloc+readname[step],'r')
    time=Data.variables['time'][:]
    datalens[step]=len(time)
    Data.close()
    
#%%
#define functions

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

print('functions defined')

#%%

#emission time
time = 0

#emission source location(lat,lon) and emission rate (g/s)

carbonsource=(450,0)

emission_rate = 21812

#directory 

loc = 'E:\Data\Co2'
writename = '\Distancedecay.nc'

#choose windspeed/direction for time in dataset

if time < datalens[0]:
    Data= Dataset(readloc+readname[0],'r')
    new_t=time
else:
    if time < datalens[1]+datalens[0]:
        Data= Dataset(readloc+readname[1],'r')
        new_t=time-datalens[0]
    else:
        if time < datalens[2]+datalens[1]+datalens[0]:
            Data= Dataset(readloc+readname[2],'r')
            new_t=time-datalens[0]-datalens[1]
        else:
            if time < datalens[3]+datalens[2]+datalens[1]+datalens[0]:
                Data= Dataset(readloc+readname[3],'r')                
                new_t=time-datalens[0]-datalens[1]-datalens[2]  
windspeed=(Data.variables['windspeed'][new_t,(900-carbonsource[0]),carbonsource[1]])
winddirection=Data.variables['winddirection'][new_t,900-carbonsource[0],carbonsource[1]]
lat= Data.variables['latitude'][:]
lon= Data.variables['longitude'][:]

#print('windspeed is ',windspeed)
#print('winddirection is ',winddirection)


#%%
#temporary change


#%%


#gaussian concentration for center point of each grid cell

grid = np.zeros((len(lat),len(lon)))
#source location(lat,lon) and emission rate (g/s)

"fix windspeed and winddirection"
#winddirection(degrees)
winddirection=90
#windspeed (m/s) 
windspeed = 4

#briggs parameters

a=0.22
b=0.20
c=0

source=[0,0]
source[0] = 900-carbonsource[0]
source[1] = carbonsource[1]

#%%

print('starting model')


for x in range(len(grid[:,0])):
    print(x)
    for y in range(len(grid[0,:])):
        #change x and y in reference to source
        xnew=(source[0]-x)*11
        ynew= (y - source[1])*5.2
        #determine distance to source
        griddirection=direction(xnew+0.000000000001,ynew+0.000000000001)
        xdist = math.cos((math.radians(griddirection-winddirection)))*(math.sqrt(xnew**2+ynew**2))
        ydist = math.sin((math.radians(griddirection-winddirection)))*(math.sqrt(xnew**2+ynew**2))
        sigma_y=a*xdist*((1+0.0001*xdist)**-0.5)
        sigma_z=b*xdist*((1+c*xdist)**-0.5)
        if sigma_z>0:
            if sigma_y>0:
                #calculate concentration at given pixel
                grid[x,y] = (emission_rate/(2*np.pi*windspeed*sigma_y*sigma_z))*math.exp(-((ydist**2)/(2*sigma_y**2)))


#%%

fig, ax = plt.subplots()
plot = plt.imshow(grid)
cbar = fig.colorbar(plot)
cbar.ax.minorticks_off()
plt.show()
#%%

"create dataset"

#create dataset
newdata= Dataset(loc+writename,'w',format="NETCDF4")
newdata= Dataset(loc+writename,'a')
#create dimensions
newlat=newdata.createDimension('latitude',len(lat))
newlon=newdata.createDimension('longitude',len(lon))
#create variables
latitudes = newdata.createVariable("latitude","f4",('latitude',))
longitudes = newdata.createVariable("longitude","f4",('longitude',))
carbons = newdata.createVariable("carbon","f8",("latitude","longitude",))
#assign data to dimensions
latitudes[:]=lat[:]
longitudes[:]=lon[:]
carbons[:]=grid[:]
#close new dataset
newdata.close()

#%%


#plt.plot(grid[450,:])
#plt.yscale("log")

